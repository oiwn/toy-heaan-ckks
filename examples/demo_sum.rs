use std::io;
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use axum::{
    Router,
    extract::State,
    http::StatusCode,
    response::Json,
    routing::{get, post},
};
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use ratatui::crossterm::{
    event::{self, Event, KeyCode},
    execute,
    terminal::{
        EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode,
        enable_raw_mode,
    },
};
use ratatui::{
    Terminal,
    backend::CrosstermBackend,
    layout::{Constraint, Direction, Layout, Rect},
    style::{Style, Stylize},
    text::{Line, Span},
    widgets::{Block, Borders, Paragraph, Wrap},
};
use reqwest;
use serde::{Deserialize, Serialize};
use toy_heaan_ckks::crypto::operations::{add_ciphertexts, decrypt, encrypt};
use toy_heaan_ckks::{
    Ciphertext, CkksEngine, EncodingParams, NaivePolyRing, Plaintext, PolyRing,
    PublicKey, SecretKey,
};

const DEGREE: usize = 8;
const SCALE_BITS: u32 = 40;
const SERVER_PORT: u16 = 3001;
const PING_INTERVAL_SECS: u64 = 2;
const CLIENT_INTERVAL_MS: u64 = 2000; // Random interval range
const MESSAGES_BEFORE_SUM: usize = 5;

#[derive(Debug)]
struct ServerState {
    ping_count: u32,
    messages: Vec<Ciphertext<NaivePolyRing<DEGREE>, DEGREE>>,
}

impl ServerState {
    fn new() -> Self {
        Self {
            ping_count: 0,
            messages: Vec::new(),
        }
    }
}

#[derive(Debug, Serialize, Deserialize)]
struct MessageRequest {
    ciphertext: SerializableCiphertext,
}

#[derive(Debug, Serialize, Deserialize)]
struct SumResponse {
    sum: SerializableCiphertext,
    count: usize,
}

#[derive(Debug, Serialize, Deserialize)]
struct ServerStatusResponse {
    message_count: usize,
    messages: Vec<String>, // Simplified view of stored messages
}

// Simple serializable wrapper for Ciphertext
#[derive(Debug, Serialize, Deserialize, Clone)]
struct SerializableCiphertext {
    c0_coeffs: Vec<u64>,
    c1_coeffs: Vec<u64>,
    scale_bits: u32,
}

#[derive(Debug, Clone)]
struct AppState {
    server_online: bool,
    last_ping_time: Option<Instant>,
    ping_response: Option<String>,
    client_messages: Vec<String>,
    client_sum_results: Vec<String>,
    server_messages: Vec<String>,
}

impl AppState {
    fn new() -> Self {
        Self {
            server_online: false,
            last_ping_time: None,
            ping_response: None,
            client_messages: Vec::new(),
            client_sum_results: Vec::new(),
            server_messages: Vec::new(),
        }
    }
}

struct CkksClientState {
    engine: CkksEngine<NaivePolyRing<DEGREE>, DEGREE>,
    secret_key: SecretKey<NaivePolyRing<DEGREE>, DEGREE>,
    public_key: PublicKey<NaivePolyRing<DEGREE>, DEGREE>,
    message_count: usize,
    encoding_params: EncodingParams<DEGREE>,
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Shared states
    let server_state = Arc::new(Mutex::new(ServerState::new()));
    let server_state_clone = server_state.clone();
    let app_state = Arc::new(Mutex::new(AppState::new()));

    // Start axum server in background thread
    let _server_handle = thread::spawn(move || {
        let rt =
            tokio::runtime::Runtime::new().expect("Failed to create tokio runtime");
        rt.block_on(async {
            start_server(server_state_clone)
                .await
                .expect("Server failed");
        });
    });

    // Give server a moment to start
    tokio::time::sleep(Duration::from_millis(100)).await;

    // Start HTTP client thread for health checks
    let app_state_health = app_state.clone();
    let _health_handle = thread::spawn(move || {
        let rt =
            tokio::runtime::Runtime::new().expect("Failed to create tokio runtime");
        rt.block_on(async {
            health_check_loop(app_state_health).await;
        });
    });

    // Start CKKS client thread
    let app_state_client = app_state.clone();
    let _client_handle = thread::spawn(move || {
        let rt =
            tokio::runtime::Runtime::new().expect("Failed to create tokio runtime");
        rt.block_on(async {
            ckks_client_loop(app_state_client).await;
        });
    });

    // Minimal, responsive TUI skeleton for CKKS demo layout
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;
    terminal.hide_cursor()?;
    terminal.clear()?;

    // Placeholders
    let mut status_note =
        String::from("Cycle 0 | Server: ❓ Unknown | Controls: [q]quit");
    let client_info: Vec<String> = vec![
        format!("Params: degree={}, scale_bits={}", DEGREE, SCALE_BITS),
        "Keys: Generated at startup".into(),
        format!("Interval: {}ms", CLIENT_INTERVAL_MS),
        format!("Sum after: {} messages", MESSAGES_BEFORE_SUM),
    ];
    let mut client_msgs: Vec<String> = vec![
        "→ ready to send values".into(),
        "→ waiting for start".into(),
    ];
    let mut store_view: Vec<String> =
        vec!["Stored: 0 items".into(), "k=... c0=[..] c1=[..]".into()];
    let logs: Vec<String> = vec![
        "[demo] TUI initialized".into(),
        format!("[server] Starting on port {}", SERVER_PORT),
        "[server] /ping endpoint ready".into(),
        "[server] /msg POST endpoint ready".into(),
        "[server] /sum GET endpoint ready".into(),
        "[server] /status GET endpoint ready".into(),
    ];

    // Tiny spinner just to show updates
    let spinner = ["|", "/", "-", "\\"];
    let mut spin_idx = 0usize;
    let mut last_tick = Instant::now();

    loop {
        // 1) Input first for instant 'q'
        if event::poll(Duration::from_millis(10))? {
            if let Event::Key(key) = event::read()? {
                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => break,
                    _ => {}
                }
            }
        }

        // 2) Lightweight periodic update for "alive" feeling + server status
        if last_tick.elapsed() >= Duration::from_millis(50) {
            spin_idx = (spin_idx + 1) % spinner.len();

            // Update server status, client messages, and server messages from app_state
            let (server_status, new_client_msgs, new_server_msgs) = {
                let app_state = app_state.lock().expect("Failed to lock app state");
                let status = if app_state.server_online {
                    if let Some(response) = &app_state.ping_response {
                        format!("✅ Online ({})", response)
                    } else {
                        "✅ Online".to_string()
                    }
                } else {
                    "❌ Offline".to_string()
                };
                (
                    status,
                    app_state.client_messages.clone(),
                    app_state.server_messages.clone(),
                )
            };

            // Update client messages display
            if !new_client_msgs.is_empty() {
                client_msgs = new_client_msgs;
            }

            // Update server messages display
            store_view = if new_server_msgs.is_empty() {
                vec!["No encrypted messages stored".into()]
            } else {
                new_server_msgs
            };

            status_note = format!(
                "Cycle 0 | Server: {} {} | Controls: [q]quit",
                server_status, spinner[spin_idx]
            );
            last_tick = Instant::now();
        }

        // 3) Draw
        terminal.draw(|f| {
            draw_ui(
                f,
                f.area(),
                &status_note,
                &client_info,
                &client_msgs,
                &store_view,
                &logs,
            )
        })?;
    }

    // Restore terminal
    terminal.show_cursor()?;
    disable_raw_mode()?;
    execute!(terminal.backend_mut(), LeaveAlternateScreen)?;

    // Note: In a real app we'd gracefully shutdown the server thread
    // For this demo, the process exit will clean up the thread

    Ok(())
}

async fn start_server(
    state: Arc<Mutex<ServerState>>,
) -> Result<(), Box<dyn std::error::Error>> {
    let app = Router::new()
        .route("/ping", get(ping_handler))
        .route("/msg", post(msg_handler))
        .route("/sum", get(sum_handler))
        .route("/status", get(status_handler))
        .with_state(state);

    let listener =
        tokio::net::TcpListener::bind(format!("127.0.0.1:{}", SERVER_PORT))
            .await
            .expect("Failed to bind server");

    axum::serve(listener, app).await.expect("Server failed");

    Ok(())
}

async fn ping_handler(
    State(state): State<Arc<Mutex<ServerState>>>,
) -> Json<String> {
    let mut state = state.lock().expect("Failed to lock server state");
    state.ping_count += 1;
    Json(format!("pong ({})", state.ping_count))
}

async fn msg_handler(
    State(state): State<Arc<Mutex<ServerState>>>,
    Json(request): Json<MessageRequest>,
) -> Result<Json<String>, StatusCode> {
    // Convert serializable format back to Ciphertext
    let ciphertext = serializable_to_ciphertext(&request.ciphertext)
        .map_err(|_| StatusCode::BAD_REQUEST)?;

    let mut state = state.lock().expect("Failed to lock server state");
    state.messages.push(ciphertext);

    Ok(Json(format!(
        "Message stored. Total: {}",
        state.messages.len()
    )))
}

async fn sum_handler(
    State(state): State<Arc<Mutex<ServerState>>>,
) -> Result<Json<SumResponse>, StatusCode> {
    let state = state.lock().expect("Failed to lock server state");

    if state.messages.is_empty() {
        return Err(StatusCode::NOT_FOUND);
    }

    // Perform homomorphic addition
    let first = &state.messages[0];
    let mut sum = Ciphertext {
        c0: first.c0.clone(),
        c1: first.c1.clone(),
        scale_bits: first.scale_bits,
    };

    for i in 1..state.messages.len() {
        sum = add_ciphertexts(&sum, &state.messages[i])
            .expect("Failed to add ciphertexts");
    }

    let serializable_sum = ciphertext_to_serializable(&sum);

    Ok(Json(SumResponse {
        sum: serializable_sum,
        count: state.messages.len(),
    }))
}

async fn status_handler(
    State(state): State<Arc<Mutex<ServerState>>>,
) -> Json<ServerStatusResponse> {
    let state = state.lock().expect("Failed to lock server state");

    // Create simplified view of stored messages for display
    let messages = state
        .messages
        .iter()
        .enumerate()
        .map(|(i, ct)| {
            format!(
                "msg_{}: scale_bits={}, c0[0]={}",
                i + 1,
                ct.scale_bits,
                ct.c0.to_coeffs()[0]
            )
        })
        .collect();

    Json(ServerStatusResponse {
        message_count: state.messages.len(),
        messages,
    })
}

// Helper functions for serialization
fn ciphertext_to_serializable(
    ct: &Ciphertext<NaivePolyRing<DEGREE>, DEGREE>,
) -> SerializableCiphertext {
    // Convert i64 to u64 for JSON serialization (simple cast for demo)
    let c0_i64 = ct.c0.to_coeffs();
    let c1_i64 = ct.c1.to_coeffs();

    SerializableCiphertext {
        c0_coeffs: c0_i64.iter().map(|&x| x as u64).collect(),
        c1_coeffs: c1_i64.iter().map(|&x| x as u64).collect(),
        scale_bits: ct.scale_bits,
    }
}

fn serializable_to_ciphertext(
    s: &SerializableCiphertext,
) -> Result<Ciphertext<NaivePolyRing<DEGREE>, DEGREE>, &'static str> {
    if s.c0_coeffs.len() != DEGREE || s.c1_coeffs.len() != DEGREE {
        return Err("Invalid coefficient array length");
    }

    // We need a context (modulus) - for now use a simple one
    // In a real app this would come from the CKKS engine
    let modulus = 741507920154517877; // Same as ckks_naive.rs

    // Convert u64 back to i64 coefficients
    let c0_coeffs: [i64; DEGREE] = s
        .c0_coeffs
        .iter()
        .map(|&x| x as i64)
        .collect::<Vec<_>>()
        .try_into()
        .map_err(|_| "Failed to convert c0 coeffs")?;

    let c1_coeffs: [i64; DEGREE] = s
        .c1_coeffs
        .iter()
        .map(|&x| x as i64)
        .collect::<Vec<_>>()
        .try_into()
        .map_err(|_| "Failed to convert c1 coeffs")?;

    let c0 = NaivePolyRing::from_coeffs(&c0_coeffs, &modulus);
    let c1 = NaivePolyRing::from_coeffs(&c1_coeffs, &modulus);

    Ok(Ciphertext {
        c0,
        c1,
        scale_bits: s.scale_bits,
    })
}

async fn health_check_loop(app_state: Arc<Mutex<AppState>>) {
    let client = reqwest::Client::new();
    let ping_url = format!("http://127.0.0.1:{}/ping", SERVER_PORT);
    let status_url = format!("http://127.0.0.1:{}/status", SERVER_PORT);

    loop {
        tokio::time::sleep(Duration::from_secs(PING_INTERVAL_SECS)).await;

        let start_time = Instant::now();

        // Check server health with ping
        let server_online = match client.get(&ping_url).send().await {
            Ok(response) => {
                if response.status().is_success() {
                    match response.text().await {
                        Ok(body) => {
                            let mut state =
                                app_state.lock().expect("Failed to lock app state");
                            state.server_online = true;
                            state.last_ping_time = Some(start_time);
                            state.ping_response =
                                Some(body.trim_matches('"').to_string());
                            true
                        }
                        Err(_) => {
                            let mut state =
                                app_state.lock().expect("Failed to lock app state");
                            state.server_online = false;
                            state.ping_response = None;
                            false
                        }
                    }
                } else {
                    let mut state =
                        app_state.lock().expect("Failed to lock app state");
                    state.server_online = false;
                    state.ping_response = None;
                    false
                }
            }
            Err(_) => {
                let mut state = app_state.lock().expect("Failed to lock app state");
                state.server_online = false;
                state.ping_response = None;
                false
            }
        };

        // If server is online, also get server status for encrypted store
        if server_online {
            match client.get(&status_url).send().await {
                Ok(response) => {
                    if let Ok(server_status) =
                        response.json::<ServerStatusResponse>().await
                    {
                        let mut state =
                            app_state.lock().expect("Failed to lock app state");
                        state.server_messages = if server_status.messages.is_empty()
                        {
                            vec![format!(
                                "📦 Storage: {} messages",
                                server_status.message_count
                            )]
                        } else {
                            let mut msgs = vec![format!(
                                "📦 Storage: {} messages",
                                server_status.message_count
                            )];
                            msgs.extend(
                                server_status
                                    .messages
                                    .iter()
                                    .map(|m| format!("🔐 {}", m)),
                            );
                            msgs
                        };
                    }
                }
                Err(_) => {
                    // Ignore status fetch errors - not critical
                }
            }
        }
    }
}

fn draw_ui(
    f: &mut ratatui::Frame,
    area: Rect,
    status_line: &str,
    client_info: &[String],
    client_msgs: &[String],
    encrypted_store: &[String],
    server_logs: &[String],
) {
    // Status (1 line), then top/bottom split
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1),
            Constraint::Percentage(70),
            Constraint::Percentage(29),
        ])
        .split(area);

    // Status line
    let status =
        Paragraph::new(Line::from(vec![Span::raw(status_line.to_string())]));
    f.render_widget(status, chunks[0]);

    // Top split: left/right
    let top_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(50), Constraint::Percentage(50)])
        .split(chunks[1]);

    // Left: Client Info (fixed height) + Client Messages (fill)
    let left_split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(7), Constraint::Min(0)])
        .split(top_chunks[0]);

    let info_block = Block::default()
        .title(Span::styled(" Client Info ", Style::default().bold()))
        .borders(Borders::ALL);
    let info_par = Paragraph::new(to_lines(client_info))
        .block(info_block)
        .wrap(Wrap { trim: true });
    f.render_widget(info_par, left_split[0]);

    let client_block = Block::default()
        .title(Span::styled(" Client Messages ", Style::default().bold()))
        .borders(Borders::ALL);
    let client_par = Paragraph::new(to_lines(client_msgs))
        .block(client_block)
        .wrap(Wrap { trim: true });
    f.render_widget(client_par, left_split[1]);

    // Right: Encrypted Store placeholder
    let enc_block = Block::default()
        .title(Span::styled(" Encrypted Store ", Style::default().bold()))
        .borders(Borders::ALL);
    let enc_par = Paragraph::new(to_lines(encrypted_store))
        .block(enc_block)
        .wrap(Wrap { trim: true });
    f.render_widget(enc_par, top_chunks[1]);

    // Bottom: Logs
    let log_block = Block::default()
        .title(Span::styled(" Logs ", Style::default().bold()))
        .borders(Borders::ALL);
    let log_par = Paragraph::new(to_lines(server_logs))
        .block(log_block)
        .wrap(Wrap { trim: true });
    f.render_widget(log_par, chunks[2]);
}

fn to_lines(items: &[String]) -> Vec<Line<'_>> {
    items.iter().map(|s| Line::from(Span::raw(s))).collect()
}

async fn ckks_client_loop(app_state: Arc<Mutex<AppState>>) {
    // Log client thread start
    {
        let mut app_state = app_state.lock().expect("Failed to lock app state");
        app_state
            .client_messages
            .insert(0, "🔧 CKKS client thread starting...".into());
    }

    // Initialize CKKS client
    let mut rng = ChaCha20Rng::from_seed([42u8; 32]); // Fixed seed for reproducibility

    // Setup CKKS engine with NaivePolyRing
    let modulus = 741507920154517877u64; // Same as ckks_naive.rs

    {
        let mut app_state = app_state.lock().expect("Failed to lock app state");
        app_state
            .client_messages
            .insert(0, "🔧 Creating CKKS engine...".into());
    }

    let engine = match CkksEngine::<NaivePolyRing<DEGREE>, DEGREE>::builder()
        .error_variance(3.2)
        .hamming_weight(DEGREE / 2) // Use same as working example
        .build_naive(modulus, SCALE_BITS)
    {
        Ok(engine) => engine,
        Err(e) => {
            let mut app_state = app_state.lock().expect("Failed to lock app state");
            app_state
                .client_messages
                .insert(0, format!("❌ CKKS engine creation failed: {}", e));
            return;
        }
    };

    {
        let mut app_state = app_state.lock().expect("Failed to lock app state");
        app_state
            .client_messages
            .insert(0, "🔑 Generating keys...".into());
    }

    // Generate keys
    let secret_key = match engine.generate_secret_key(&mut rng) {
        Ok(sk) => sk,
        Err(e) => {
            let mut app_state = app_state.lock().expect("Failed to lock app state");
            app_state
                .client_messages
                .insert(0, format!("❌ Secret key generation failed: {}", e));
            return;
        }
    };

    let public_key = match engine.generate_public_key(&secret_key, &mut rng) {
        Ok(pk) => pk,
        Err(e) => {
            let mut app_state = app_state.lock().expect("Failed to lock app state");
            app_state
                .client_messages
                .insert(0, format!("❌ Public key generation failed: {}", e));
            return;
        }
    };

    // Create encoding params
    let encoding_params = match EncodingParams::new(1 << SCALE_BITS) {
        Ok(params) => params,
        Err(e) => {
            let mut app_state = app_state.lock().expect("Failed to lock app state");
            app_state
                .client_messages
                .insert(0, format!("❌ Encoding params creation failed: {}", e));
            return;
        }
    };

    let mut client_state = CkksClientState {
        engine,
        secret_key,
        public_key,
        message_count: 0,
        encoding_params,
    };

    let client = reqwest::Client::new();
    let msg_url = format!("http://127.0.0.1:{}/msg", SERVER_PORT);
    let sum_url = format!("http://127.0.0.1:{}/sum", SERVER_PORT);

    {
        let mut app_state = app_state.lock().expect("Failed to lock app state");
        app_state
            .client_messages
            .insert(0, "✅ CKKS client initialized successfully!".into());
    }

    // Give server time to start
    tokio::time::sleep(Duration::from_millis(500)).await;

    {
        let mut app_state = app_state.lock().expect("Failed to lock app state");
        app_state
            .client_messages
            .insert(0, "🚀 Starting message sending loop...".into());
    }

    loop {
        // Wait random interval
        let wait_time = CLIENT_INTERVAL_MS + (rng.random::<u64>() % 1000);
        tokio::time::sleep(Duration::from_millis(wait_time)).await;

        // Generate random value in [-10, 10]
        let random_value = rng.random_range(-10.0..=10.0);

        // Log attempt
        {
            let mut app_state = app_state.lock().expect("Failed to lock app state");
            app_state
                .client_messages
                .insert(0, format!("📝 Encrypting value: {:.2}...", random_value));
        }

        // Encrypt the value
        match encrypt_and_send_value(&client_state, &client, &msg_url, random_value)
            .await
        {
            Ok(response) => {
                client_state.message_count += 1;

                // Update app state with new message
                {
                    let mut app_state =
                        app_state.lock().expect("Failed to lock app state");
                    app_state.client_messages.insert(
                        0,
                        format!(
                            "📤 Sent: {:.2} (encrypted) - {}",
                            random_value, response
                        ),
                    );

                    // Keep only last 10 messages
                    if app_state.client_messages.len() > 10 {
                        app_state.client_messages.truncate(10);
                    }
                }

                // Request sum after MESSAGES_BEFORE_SUM messages
                if client_state.message_count >= MESSAGES_BEFORE_SUM {
                    match request_and_decrypt_sum(&client_state, &client, &sum_url)
                        .await
                    {
                        Ok(sum_result) => {
                            let mut app_state =
                                app_state.lock().expect("Failed to lock app state");
                            app_state.client_messages.insert(
                                0,
                                format!(
                                    "🧮 Sum of {} values: {:.2}",
                                    MESSAGES_BEFORE_SUM, sum_result
                                ),
                            );
                            client_state.message_count = 0; // Reset counter
                        }
                        Err(e) => {
                            let mut app_state =
                                app_state.lock().expect("Failed to lock app state");
                            app_state
                                .client_messages
                                .insert(0, format!("❌ Sum request failed: {}", e));
                        }
                    }
                }
            }
            Err(e) => {
                let mut app_state =
                    app_state.lock().expect("Failed to lock app state");
                app_state
                    .client_messages
                    .insert(0, format!("❌ Send failed: {}", e));
            }
        }
    }
}

async fn encrypt_and_send_value(
    client_state: &CkksClientState,
    client: &reqwest::Client,
    url: &str,
    value: f64,
) -> Result<String, Box<dyn std::error::Error>> {
    // Encode and encrypt the single value
    let values = [value, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; // Pad to DEGREE length
    let encoded_coeffs =
        toy_heaan_ckks::encode(&values, &client_state.encoding_params)
            .expect("Failed to encode");

    // Create polynomial from encoded coefficients
    let modulus = 741507920154517877u64;
    let poly = NaivePolyRing::from_coeffs(&encoded_coeffs, &modulus);
    let plaintext = Plaintext {
        poly,
        scale_bits: SCALE_BITS,
    };

    let mut rng = ChaCha20Rng::from_seed([42u8; 32]);
    let ciphertext = encrypt(
        &plaintext,
        &client_state.public_key,
        &mut rng,
        &modulus,
        3.2,
    )
    .expect("Failed to encrypt");

    // Convert to serializable format
    let serializable_ct = ciphertext_to_serializable(&ciphertext);
    let request = MessageRequest {
        ciphertext: serializable_ct,
    };

    // Send POST request
    let response = client.post(url).json(&request).send().await?.text().await?;

    Ok(response.trim_matches('"').to_string())
}

async fn request_and_decrypt_sum(
    client_state: &CkksClientState,
    client: &reqwest::Client,
    url: &str,
) -> Result<f64, Box<dyn std::error::Error>> {
    // Request sum from server
    let response: SumResponse = client.get(url).send().await?.json().await?;

    // Convert back to ciphertext and decrypt
    let sum_ciphertext = serializable_to_ciphertext(&response.sum)
        .map_err(|e| format!("Failed to deserialize sum: {}", e))?;

    let decrypted_plaintext = decrypt(&sum_ciphertext, &client_state.secret_key);

    // Decode back to float
    let coeffs = decrypted_plaintext.poly.to_coeffs();
    let decoded_values =
        toy_heaan_ckks::decode(&coeffs, &client_state.encoding_params)
            .expect("Failed to decode sum");

    Ok(decoded_values[0]) // Return first element
}
