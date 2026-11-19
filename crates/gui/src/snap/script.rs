use eframe::egui::{self, Key, Modifiers, PointerButton, Pos2, Rect};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Crop {
    #[default]
    Full,
    Viewport,
}

impl Crop {
    pub fn parse(value: &str) -> Result<Self, String> {
        match value {
            "full" => Ok(Self::Full),
            "viewport" => Ok(Self::Viewport),
            _ => Err(format!("expected full or viewport, got {value:?}")),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Position {
    Window(f32, f32),
    Viewport(f32, f32),
}

impl Position {
    fn parse(value: &str) -> Result<Self, String> {
        let (fraction, pair) = value
            .strip_prefix("vp:")
            .map_or((false, value), |pair| (true, pair));
        let (x, y) = pair
            .split_once(',')
            .ok_or_else(|| format!("expected x,y or vp:fx,fy, got {value:?}"))?;
        let (x, y) = (number(x)?, number(y)?);
        if fraction {
            if !(0.0..=1.0).contains(&x) || !(0.0..=1.0).contains(&y) {
                return Err("viewport fractions must be between 0 and 1".into());
            }
            Ok(Self::Viewport(x, y))
        } else {
            Ok(Self::Window(x, y))
        }
    }

    pub(super) fn resolve(self, viewport: Option<Rect>, scale: f32) -> Result<Pos2, String> {
        match self {
            Self::Window(x, y) => Ok(Pos2::new(x / scale, y / scale)),
            Self::Viewport(x, y) => {
                let r = viewport.ok_or("no 3D viewport on screen")?;
                Ok(Pos2::new(
                    r.left() + x * r.width(),
                    r.top() + y * r.height(),
                ))
            }
        }
    }
}

#[derive(Debug)]
pub(super) enum Step {
    Settle(f32),
    Wait(u32),
    Module(crate::Module),
    Snap(String, Option<Crop>),
    Key(Key, Modifiers),
    Text(String),
    Click(Position, PointerButton),
    Label(String),
    Hover(Position),
    Drag(Modifiers, PointerButton, Position, egui::Vec2, u32),
    Scroll(f32, Option<Position>),
    Pinch(f32, Option<Position>),
    Trace(String),
    End,
}

pub(super) struct Instruction {
    pub source: String,
    pub step: Step,
}

pub(super) fn parse(input: &str) -> Result<Vec<Instruction>, String> {
    let mut script = Vec::new();
    for (index, line) in input.lines().enumerate() {
        let at_line = |error: String| format!("line {}: {error}", index + 1);
        for tokens in tokenize(line).map_err(at_line)? {
            let source = tokens
                .iter()
                .map(|token| {
                    if token.quoted {
                        format!("{:?}", token.value)
                    } else {
                        token.value.clone()
                    }
                })
                .collect::<Vec<_>>()
                .join(" ");
            let step = parse_step(&tokens).map_err(at_line)?;
            script.push(Instruction { source, step });
        }
    }
    Ok(script)
}

#[derive(Default)]
struct Token {
    value: String,
    quoted: bool,
}

fn tokenize(line: &str) -> Result<Vec<Vec<Token>>, String> {
    let mut commands = Vec::new();
    let mut tokens = Vec::new();
    let mut token = Token::default();
    let mut chars = line.chars();
    let mut in_quote = false;
    while let Some(ch) = chars.next() {
        if in_quote {
            match ch {
                '"' => in_quote = false,
                '\\' => token.value.push(chars.next().ok_or("unfinished escape")?),
                _ => token.value.push(ch),
            }
        } else {
            match ch {
                '"' if token.value.is_empty() && !token.quoted => {
                    token.quoted = true;
                    in_quote = true;
                }
                '"' => return Err("put quotes around a whole argument".into()),
                '#' => break,
                ';' => {
                    push_token(&mut tokens, &mut token);
                    if !tokens.is_empty() {
                        commands.push(std::mem::take(&mut tokens));
                    }
                }
                ch if ch.is_whitespace() => push_token(&mut tokens, &mut token),
                _ if token.quoted => return Err("expected a space after a quoted argument".into()),
                _ => token.value.push(ch),
            }
        }
    }
    if in_quote {
        return Err("unterminated quote".into());
    }
    push_token(&mut tokens, &mut token);
    if !tokens.is_empty() {
        commands.push(tokens);
    }
    Ok(commands)
}

fn push_token(tokens: &mut Vec<Token>, token: &mut Token) {
    if !token.value.is_empty() || token.quoted {
        tokens.push(std::mem::take(token));
    }
}

fn number(value: &str) -> Result<f32, String> {
    value
        .parse()
        .ok()
        .filter(|n: &f32| n.is_finite())
        .ok_or_else(|| format!("expected a finite number, got {value:?}"))
}

fn positive(value: &str) -> Result<f32, String> {
    let parsed = number(value)?;
    if parsed > 0.0 {
        Ok(parsed)
    } else {
        Err(format!("expected a positive number, got {value:?}"))
    }
}

fn frames(value: &str) -> Result<u32, String> {
    value
        .parse()
        .map_err(|_| format!("expected a frame count, got {value:?}"))
}

fn name(value: &str) -> Result<String, String> {
    if value.is_empty()
        || !value
            .bytes()
            .all(|ch| ch.is_ascii_alphanumeric() || ch == b'-' || ch == b'_')
    {
        return Err(format!(
            "invalid capture name {value:?}, use letters, digits, hyphens or underscores"
        ));
    }
    Ok(value.into())
}

fn modified(value: &str) -> Result<(Modifiers, &str), String> {
    let mut pieces = value.split('+').peekable();
    let mut mods = Modifiers::NONE;
    while let Some(part) = pieces.next() {
        if pieces.peek().is_none() {
            return Ok((mods, part));
        }
        match part.to_ascii_lowercase().as_str() {
            "ctrl" => {
                mods.ctrl = true;
                mods.command = true;
            }
            "shift" => mods.shift = true,
            "alt" => mods.alt = true,
            _ => return Err(format!("unknown modifier {part:?}")),
        }
    }
    Err("expected a key or mouse button".into())
}

fn button(value: &str) -> Result<PointerButton, String> {
    match value.to_ascii_lowercase().as_str() {
        "lmb" => Ok(PointerButton::Primary),
        "mmb" => Ok(PointerButton::Middle),
        "rmb" => Ok(PointerButton::Secondary),
        _ => Err(format!("expected lmb, mmb or rmb, got {value:?}")),
    }
}

fn parse_step(tokens: &[Token]) -> Result<Step, String> {
    let args: Vec<_> = tokens.iter().map(|token| token.value.as_str()).collect();
    let step = match args.as_slice() {
        ["settle"] => Step::Settle(120.0),
        ["settle", seconds] => Step::Settle(positive(seconds)?.min(120.0)),
        ["wait", count] => Step::Wait(frames(count)?),
        ["module", module] => Step::Module(super::module(module)?),
        ["snap", stem] => Step::Snap(name(stem)?, None),
        ["snap", stem, crop] => Step::Snap(name(stem)?, Some(Crop::parse(crop)?)),
        ["key", chord] => {
            let (mods, key) = modified(chord)?;
            let key = Key::from_name(key)
                .or_else(|| Key::from_name(&key.to_ascii_uppercase()))
                .ok_or_else(|| format!("unknown key {key:?}"))?;
            Step::Key(key, mods)
        }
        ["text", value] => Step::Text((*value).into()),
        ["click", label] if tokens[1].quoted => Step::Label((*label).into()),
        ["click", pos] => Step::Click(Position::parse(pos)?, PointerButton::Primary),
        ["click", pos, which] => Step::Click(Position::parse(pos)?, button(which)?),
        ["hover", pos] => Step::Hover(Position::parse(pos)?),
        ["drag", chord, from, dx, dy, rest @ ..] if rest.len() <= 1 => {
            let (mods, which) = modified(chord)?;
            let count = rest.first().map_or(Ok(10), |count| frames(count))?;
            if count == 0 {
                return Err("drag needs at least one motion frame".into());
            }
            Step::Drag(
                mods,
                button(which)?,
                Position::parse(from)?,
                egui::vec2(number(dx)?, number(dy)?),
                count,
            )
        }
        ["scroll", dy, pos @ ..] if pos.len() <= 1 => Step::Scroll(
            number(dy)?,
            pos.first().map(|p| Position::parse(p)).transpose()?,
        ),
        ["pinch", factor, pos @ ..] if pos.len() <= 1 => Step::Pinch(
            positive(factor)?,
            pos.first().map(|p| Position::parse(p)).transpose()?,
        ),
        ["trace", stem] => Step::Trace(name(stem)?),
        ["end"] => Step::End,
        _ => return Err(format!("invalid snap step {:?}", args.join(" "))),
    };
    Ok(step)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_scripts_and_reports_bad_steps() {
        let script = parse("settle 8; click \"Label #1; A\" # ignore\n drag ctrl+shift+mmb vp:0.5,0.4 20 -5 4\n key shift+1").unwrap();
        assert!(matches!(script[0].step, Step::Settle(8.0)));
        assert!(matches!(&script[1].step, Step::Label(label) if label == "Label #1; A"));
        assert!(
            matches!(script[2].step, Step::Drag(mods, PointerButton::Middle, _, _, 4) if mods.ctrl && mods.shift)
        );
        assert!(matches!(script[3].step, Step::Key(Key::Num1, mods) if mods.shift));
        for invalid in [
            "click \"missing",
            "drag lmb 0,0 1 2 0",
            "snap ../x",
            "pinch 0",
            "key no-such-key",
            "hover vp:2,0",
        ] {
            assert!(parse(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn resolves_pixels_and_viewport_fractions() {
        let rect = Rect::from_min_max(Pos2::new(100.0, 20.0), Pos2::new(300.0, 120.0));
        assert_eq!(
            Position::parse("vp:0.25,0.5")
                .unwrap()
                .resolve(Some(rect), 2.0)
                .unwrap(),
            Pos2::new(150.0, 70.0)
        );
        assert_eq!(
            Position::parse("200,120")
                .unwrap()
                .resolve(None, 2.0)
                .unwrap(),
            Pos2::new(100.0, 60.0)
        );
        assert!(Position::parse("vp:0,0")
            .unwrap()
            .resolve(None, 1.0)
            .is_err());
    }
}
