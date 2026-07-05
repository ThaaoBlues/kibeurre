/*
All credits go to Gemini for the TUI implementation and design.
*/

use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{
    backend::{Backend, CrosstermBackend},
    layout::{Constraint, Direction, Layout},
    style::{Color, Modifier, Style},
    widgets::{Block, Borders, Paragraph, Wrap},
    Terminal,
};
use std::io;

// Pulling in your modules and constants
use crate::parameters::{d_u, d_v, eta_1, eta_2, k, m, n, q};


use crate::core::{encrypt, decrypt, generate_seed_vector, generate_A_from_seed, generate_noise_polyvector, compute_t, EncryptedMessage};
use crate::format_utils::{string_to_vectors, vectors_to_string};


/// Tracks the lifecycle events of the encryption run
struct CryptoStep {
    name: String,
    details: String,
}

struct AppState {
    input_buffer: String,
    steps_log: Vec<CryptoStep>,
    decrypted_result: String,
    success: Option<bool>,
    logs_scroll: u16,
}

impl AppState {
    fn new() -> Self {
        AppState {
            input_buffer: String::new(),
            steps_log: Vec::new(),
            decrypted_result: String::new(),
            success: None,
            logs_scroll: 0,
        }
    }

    fn run_kyber_pipeline(&mut self) {
        if self.input_buffer.is_empty() {
            return;
        }
        self.steps_log.clear();
        self.logs_scroll = 0;

        // 1. Setup Phase
        self.steps_log.push(CryptoStep {
            name: "[1] Key Generation".to_string(),
            details: format!("Generated seed vector size {n}.\nMatrix A constructed via Shake128 (dimension: {k}x{k})."),
        });
        
        let seed_vector = generate_seed_vector();
        let A = generate_A_from_seed(seed_vector);
        let mut s = generate_noise_polyvector(eta_1);
        let e = generate_noise_polyvector(eta_2);
        let t = compute_t(A.clone(), s.clone(), e);

        
        let msg = string_to_vectors(&self.input_buffer);

        self.steps_log.push(CryptoStep {
            name: "[1.1] Message Polynomial Mapping".to_string(),
            details: format!("Input string as vectors of size {n} generated: {:?}", msg),
        });

        
        self.steps_log.push(CryptoStep {
            name: "[2] Message Polynomial Mapping".to_string(),
            details: format!("Input string processed and mapped onto an R_q polynomial degree vector."),
        });

        // 3. Encryption
        self.steps_log.push(CryptoStep {
            name: "[3] Encryption (Capsule)".to_string(),
            details: format!("Generated error vectors (η1={eta_1}, η2={eta_2}).\nComputed Vector pair (u, v)."),
        });

        let mut encrypted_chunks = Vec::new();
        for (i, chunk) in msg.iter().enumerate() {
            let r = generate_noise_polyvector(eta_1);
            encrypted_chunks.push(encrypt(A, t, r, chunk.clone()));
        }

        /*self.steps_log.push(CryptoStep {
            name: "[3.1] Encryption Output".to_string(),
            details: format!("Encrypted message chunks generated: {:?}", encrypted_chunks),
        });*/


        // 4. Decryption
        self.steps_log.push(CryptoStep {
            name: "[4] Decryption".to_string(),
            details: "Running error correction round operation: m' = round(v - s^T * u)".to_string(),
        });


        let mut decrypted_chunks = Vec::new();

        for encrypted in &mut encrypted_chunks {
            let decrypted_msg = decrypt(encrypted, &mut s);
            decrypted_chunks.push(decrypted_msg);
        }

        self.steps_log.push(CryptoStep {
            name: "[4.1] Decryption Output".to_string(),
            details: format!("Decrypted message chunks generated: {:?}", decrypted_chunks),
        });

        let decrypted_msg = vectors_to_string(decrypted_chunks);


        self.success = Some(decrypted_msg == self.input_buffer);
        self.decrypted_result = decrypted_msg.clone();

        self.steps_log.push(CryptoStep {
            name: "[5] Verification".to_string(),
            details: format!("Decrypted message: \"{}\"\nOriginal message: \"{}\"\nMatch: {}", decrypted_msg, self.input_buffer, self.success.unwrap()),
        });

    }
}

pub fn mainloop() -> Result<(), Box<dyn std::error::Error>> {
    // Setup Terminal
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    let mut app_state = AppState::new();
    let res = run_app(&mut terminal, &mut app_state);

    // Restore Terminal
    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    if let Err(err) = res {
        println!("{:?}", err);
    }
    Ok(())
}

fn run_app<B: Backend>(terminal: &mut Terminal<B>, app: &mut AppState) -> io::Result<()> {
    loop {
        terminal.draw(|f| ui(f, app));

        if let Event::Key(key) = event::read()? {
            match key.code {
                KeyCode::Esc => return Ok(()),
                KeyCode::Enter => {
                    app.run_kyber_pipeline();
                }
                KeyCode::Char(c) => {
                    app.input_buffer.push(c);
                }
                KeyCode::Backspace => {
                    app.input_buffer.pop();
                }

                // Log scrolling controls
                KeyCode::Up => {
                    if app.logs_scroll > 0 {
                        app.logs_scroll -= 1;
                    }
                }
                KeyCode::Down => {
                    // Prevent infinite scrolling downwards
                    if !app.steps_log.is_empty() && app.logs_scroll < 100 { 
                        app.logs_scroll += 1;
                    }
                }
                _ => {}
            }
        }
    }
}

fn ui(f: &mut ratatui::Frame, app: &AppState) {
    // Layout Splitting
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3), // Header title banner
            Constraint::Min(10),   // Workspace body
            Constraint::Length(3), // Bottom status bar
        ])
        .split(f.size());

    // Title Widget
    let header = Paragraph::new("KIBEURREEEEEEEE TEST TUI")
        .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD))
        .block(Block::default().borders(Borders::ALL).border_style(Style::default().fg(Color::Cyan)));
    f.render_widget(header, chunks[0]);

    // Workspace Division (Left Panel for Logs, Right Panel for Config Parameters)
    let workspace_chunks = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(70), Constraint::Percentage(30)])
        .split(chunks[1]);

    // Left Workspace: Sub-split into input box, sequence execution steps, and output verification
    let left_chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3), // Input field
            Constraint::Min(5),    // Cryptographic pipeline traces
            Constraint::Length(4), // Decryption output window
        ])
        .split(workspace_chunks[0]);

    // Input Message Box
    let input_widget = Paragraph::new(app.input_buffer.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Enter Plaintext Message to Encrypt "));
    f.render_widget(input_widget, left_chunks[0]);

    // Render Steps Execution Stream
    let mut logs_text = String::new();
    if app.steps_log.is_empty() {
        logs_text = "\n Press [Enter] to run the message through your Kyber engine implementation pipeline.".to_string();
    } else {
        for step in &app.steps_log {
            logs_text.push_str(&format!("● {}\n  {}\n\n", step.name, step.details));
        }
    }


    // Output Result Window
    let status_color = match app.success {
        Some(true) => Color::Green,
        Some(false) => Color::Red,
        None => Color::DarkGray,
    };
    let result_widget = Paragraph::new(app.decrypted_result.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Decrypted Output Validation ").border_style(Style::default().fg(status_color)));
    f.render_widget(result_widget, left_chunks[2]);

    // Right Workspace: Parameter Inspection Side Panel
    let param_text = format!(
        "\n \
         • Degree (n): {}\n \
         • Modules (k): {}\n \
         • Field (q): {}\n \
         • η1 (Noise): {}\n \
         • η2 (Noise): {}\n \
         • d_u (Bits): {}\n \
         • d_v (Bits): {}\n \
         • Ring Mod (m): {}\n\n \
         ",
        n, k, q, eta_1, eta_2, d_u, d_v, m
    );
    let params_widget = Paragraph::new(param_text)
        .block(Block::default().borders(Borders::ALL).title(" Kyber Configuration ").border_style(Style::default().fg(Color::Magenta)))
        .wrap(Wrap { trim: true });
    f.render_widget(params_widget, workspace_chunks[1]);

    // Bottom Quick Help Bar
    let help_widget = Paragraph::new(" Type message text | [Enter] Simulate Cryptosystem Run | [Esc] Exit Interactive Mode")
        .style(Style::default().fg(Color::DarkGray));
    f.render_widget(help_widget, chunks[2]);


    // Render Steps Execution Stream
    let mut logs_text = String::new();
    if app.steps_log.is_empty() {
        logs_text = "\n Press [Enter] to run the message through your Kyber engine implementation pipeline.".to_string();
    } else {
        for step in &app.steps_log {
            logs_text.push_str(&format!("● {}\n  {}\n\n", step.name, step.details));
        }
    }
    
    // Add scroll tuple parameter: (row_scroll, col_scroll)
    let logs_widget = Paragraph::new(logs_text)
        .block(Block::default()
            .borders(Borders::ALL)
            .title(format!(" Execution Steps & Matrix Operations [Scroll: {}] ", app.logs_scroll))
        )
        .wrap(Wrap { trim: true })
        .scroll((app.logs_scroll, 0)); // <--- Apply the vertical scroll state modifier here
        
    f.render_widget(logs_widget, left_chunks[1]);
}