

const N: i32 = 1;
const runtime: f64 = 1.0; // s
const dt: f64 = 1e-6; // s

const radius: f64 = 1e-9; // m
const cell_viscosity: f64 = 7e-4; // kg / s*m
const gamma: f64 = radius * 6.0 * std::f64::consts::PI * cell_viscosity; // kg/s

const kb: f64 = 1.38e-23; // J/T
const T: f64 = 310.15; // K

fn simulate_particle(init_position: Vec<f64>) -> Vec<Vec<f64>> {
    vec![vec![0.0, 0.0, 0.0]]
}









fn main() {
    let mut trajectories = Vec::new();
    for _ in 0..N {
        let mut init_position = vec![0.0, 0.0, 0.0];
        let trajectory: Vec<Vec<f64>> = simulate_particle(init_position);
        trajectories.push(trajectory);
    }
    println!("trajectories: {:?}\n", trajectories);
}
