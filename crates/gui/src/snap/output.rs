use super::{Capture, SnapError};
use font8x8::{UnicodeFonts, BASIC_FONTS};
use image::{codecs::gif, imageops, Delay, Frame, Rgba, RgbaImage};
use std::{fs::File, path::Path};

/// A contact sheet uses evenly spaced frames, including the first and last.
pub(super) fn sheet(
    dir: &Path,
    frames: &[&Capture],
    file: &str,
    columns: u32,
    limit: usize,
) -> Result<(), SnapError> {
    let count = frames.len().min(limit);
    if count == 0 {
        return Err(SnapError::Input(format!("no frames to write to {file}")));
    }
    const CELL_WIDTH: u32 = 352;
    const CELL_HEIGHT: u32 = 244;
    let rows = (count as u32).div_ceil(columns);
    let mut image = RgbaImage::from_pixel(
        columns * CELL_WIDTH,
        rows * CELL_HEIGHT,
        Rgba([22, 25, 31, 255]),
    );
    for i in 0..count {
        let index = if count == 1 {
            0
        } else {
            i * (frames.len() - 1) / (count - 1)
        };
        let capture = frames[index];
        let frame = read_frame(dir, capture)?.to_rgba8();
        let thumb = imageops::thumbnail(&frame, CELL_WIDTH - 8, CELL_HEIGHT - 22);
        let x = (i as u32 % columns) * CELL_WIDTH;
        let y = (i as u32 / columns) * CELL_HEIGHT;
        imageops::overlay(
            &mut image,
            &thumb,
            i64::from(x + (CELL_WIDTH - thumb.width()) / 2),
            i64::from(y + 4 + (CELL_HEIGHT - 22 - thumb.height()) / 2),
        );
        let label = if let Some(view) = &capture.viewport {
            format!(
                "{} yaw {:.1} pitch {:.1} dist {:.2}",
                capture.frame,
                view.yaw.to_degrees(),
                view.pitch.to_degrees(),
                view.distance
            )
        } else {
            format!("frame {}", capture.frame)
        };
        caption(&mut image, x + 6, y + CELL_HEIGHT - 13, &label);
    }
    image.save(dir.join(file)).map_err(SnapError::output)
}

fn caption(image: &mut RgbaImage, x: u32, y: u32, text: &str) {
    for (index, ch) in text.chars().enumerate() {
        if let Some(glyph) = BASIC_FONTS.get(ch) {
            for (row, bits) in glyph.into_iter().enumerate() {
                for col in 0..8 {
                    if bits & (1 << col) != 0 {
                        let pos = x + index as u32 * 8 + col;
                        if pos < image.width() {
                            image.put_pixel(pos, y + row as u32, Rgba([236, 240, 246, 255]));
                        }
                    }
                }
            }
        }
    }
}

fn read_frame(dir: &Path, capture: &Capture) -> Result<image::DynamicImage, SnapError> {
    image::open(dir.join(&capture.file)).map_err(SnapError::output)
}

pub(super) fn animation(dir: &Path, frames: &[&Capture], file: &str) -> Result<(), SnapError> {
    let mut encoder =
        gif::GifEncoder::new(File::create(dir.join(file)).map_err(SnapError::output)?);
    encoder
        .set_repeat(gif::Repeat::Infinite)
        .map_err(SnapError::output)?;
    for frame in frames {
        let image = read_frame(dir, frame)?.to_rgba8();
        encoder
            .encode_frame(Frame::from_parts(
                image,
                0,
                0,
                Delay::from_numer_denom_ms(1000, 60),
            ))
            .map_err(SnapError::output)?;
    }
    Ok(())
}
