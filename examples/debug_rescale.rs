use toy_heaan_ckks::{NaivePolyRing, PolyRing};

const DEGREE: usize = 8;
const MODULUS: u64 = (1u64 << 50) - 27; // Large modulus for headroom

fn main() {
    println!("🔬 Testing modular rescaling fix");
    println!("   Modulus: {}", MODULUS);
    println!("   Half modulus: {}", MODULUS / 2);

    // Test our problematic coefficient from mul_try.rs
    let raw_coeff = 1125899906841473u64;
    println!("\n🔍 Analyzing coefficient: {}", raw_coeff);

    // Check if this represents a negative value
    let half_modulus = MODULUS / 2;
    if raw_coeff > half_modulus {
        let negative_value = -((MODULUS - raw_coeff) as i64);
        println!("   This represents: {} (negative)", negative_value);
        println!(
            "   After scaling by 1024: {:.6}",
            negative_value as f64 / 1024.0
        );
    } else {
        println!("   This represents: {} (positive)", raw_coeff);
        println!("   After scaling by 1024: {:.6}", raw_coeff as f64 / 1024.0);
    }

    // Test the to_coeffs conversion
    let mut poly = NaivePolyRing::<DEGREE>::with_modulus(MODULUS);
    poly.coeffs[0] = raw_coeff;

    let signed_coeffs = poly.to_coeffs();
    println!("   to_coeffs() gives: {}", signed_coeffs[0]);
    println!("   After scaling: {:.6}", signed_coeffs[0] as f64 / 1024.0);

    // Test some other coefficients from our debug output
    let coeffs = [2493u64, 586u64, 1125899906841473u64, 1869u64];
    println!("\n📊 All coefficients from debug:");
    for (i, &coeff) in coeffs.iter().enumerate() {
        if coeff > half_modulus {
            let negative_value = -((MODULUS - coeff) as i64);
            println!(
                "   coeffs[{}]: {} → {} → {:.6} (after scaling)",
                i,
                coeff,
                negative_value,
                negative_value as f64 / 1024.0
            );
        } else {
            println!(
                "   coeffs[{}]: {} → {:.6} (after scaling)",
                i,
                coeff,
                coeff as f64 / 1024.0
            );
        }
    }

    println!("\n✅ Expected range for results [2, 2, 3, 8]:");
    println!("   After encoding with scale 1024: [2048, 2048, 3072, 8192]");
    println!("   So these coefficients should be close to those values!");
}
