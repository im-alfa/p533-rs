use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestConfig {
    pub general: GeneralConfig,
    pub scenarios: Vec<TestScenario>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneralConfig {
    pub tolerance_default: f64,
    pub tolerance_strict: f64,
    pub year: i32,
    pub month: u8,
    pub hour: u8,
    pub ssn_normal: i32,
    pub ssn_low: i32,
    pub ssn_high: i32,
    pub frequency_default: f64,
    pub power: f64,
    pub tx_gain: f64,
    pub rx_gain: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestScenario {
    pub name: String,
    pub description: String,
    pub tx_lat: f64,
    pub tx_lon: f64,
    pub rx_lat: f64,
    pub rx_lon: f64,
    pub frequency: f64,
    pub hour: u8,
    pub ssn: i32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestResult {
    pub name: String,
    pub distance: f64,
    pub muf50: f64,
    pub bcr: f64,
    pub snr: f64,
    pub e_layer_screening_frequency: f64,
    pub median_available_receiver_power: f64,
}

pub fn load_test_config() -> TestConfig {
    let config_path = "tests/test_config.json";
    let config_content = fs::read_to_string(config_path).expect("Failed to read test config file");

    serde_json::from_str(&config_content).expect("Failed to parse test config")
}

pub fn create_test_scenarios() -> Vec<TestScenario> {
    let config = load_test_config();
    config.scenarios
}

pub fn run_p533_prediction(scenario: &TestScenario) -> TestResult {
    let config = load_test_config();

    let path = p533::create_path_prediction(
        scenario.tx_lat,
        scenario.tx_lon,
        scenario.rx_lat,
        scenario.rx_lon,
        config.general.year,
        config.general.month as i32,
        scenario.hour as i32,
        scenario.ssn,
        scenario.frequency,
        config.general.power,
        config.general.tx_gain,
        config.general.rx_gain,
    );

    TestResult {
        name: scenario.name.clone(),
        distance: path.distance,
        muf50: path.muf50,
        bcr: path.bcr,
        snr: path.snr,
        e_layer_screening_frequency: path.es,
        median_available_receiver_power: path.pr,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let scenarios = create_test_scenarios();
    let mut results = Vec::new();

    println!("running {} scenarios", scenarios.len());

    for scenario in &scenarios {
        let result = run_p533_prediction(scenario);
        println!(
            "{}: d={:.1}km, muf={:.2}MHz, bcr={:.1}%",
            result.name, result.distance, result.muf50, result.bcr
        );
        results.push(result);
    }

    // save results
    let json_output = serde_json::to_string_pretty(&results)?;
    std::fs::write("tests/rust_results.json", json_output)?;

    println!("saved to tests/rust_results.json");

    Ok(())
}
