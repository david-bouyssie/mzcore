use anyhow::*;
//use itertools::Itertools;
//use std::str::Chars;
//use std::collections::HashMap;

//use crate::chemistry::composition::*;
use crate::chemistry::table::AminoAcidTable;
use crate::chemistry::constants::{WATER_MONO_MASS, WATER_AVERAGE_MASS};

pub fn calc_aa_seq_mass(aa_seq: &str, aa_table: &AminoAcidTable, mono_mass: bool) -> Result<f64> {
    let vec = aa_seq.chars().collect::<Vec<_>>();
    let char_iter = vec.iter();
    calc_aa_seq_mass_from_chars(char_iter, aa_table, mono_mass)
}

pub fn calc_aa_seq_mass_from_chars(aa_as_chars: std::slice::Iter<char>, aa_table: &AminoAcidTable, mono_mass: bool) -> Result<f64> {

    let aa_by_code1 = &aa_table.aa_by_code1;

    let mut seq_mass = 0.0;

    for aa_as_char in aa_as_chars {
        //let aa_as_char = seq_as_bytes[char_idx] as char;
        let aa_opt = aa_by_code1.get(&aa_as_char);
        let aa = aa_opt.context(format!(
            "can't find amino acid '{}' in the provided table",aa_as_char
        ))?;

        seq_mass += if mono_mass {aa.mono_mass} else {aa.average_mass};
    }

    seq_mass += if mono_mass {WATER_MONO_MASS} else {WATER_AVERAGE_MASS};

    Ok(seq_mass)
}

/*pub fn calc_aa_seq_mass(aa_seq: &str, aa_table: &AminoAcidTable, mono_mass: bool) -> Result<f64> {

    let aa_comp = parse_aa_composition(aa_seq)?;

    let get_aa_mass = |aa_code1: char| -> Result<f64> {
        let aa = aa_table.aa_by_code1.get(&aa_code1).ok_or_else(
            || anyhow!("can't find amino acid '{}' in the provided table",aa_code1)
        )?;
        let m = if mono_mass {aa.mono_mass} else {aa.average_mass};
        Ok(m)
    };

    let seq_mass = _calc_mass(aa_comp, get_aa_mass)?;

    if mono_mass {
        Ok(seq_mass + WATER_MONO_MASS)
    } else {
        Ok(seq_mass + WATER_AVERAGE_MASS)
    }
}

fn _calc_mass<T,F>(abundance_map: HashMap<T, f32>, get_entity_mass: F) -> Result<f64> where F: Fn(T) -> Result<f64> {

    let mut mass: f64 = 0.0;
    for (entity, entity_ab) in abundance_map {
        let entity_mass = get_entity_mass(entity)?;
        mass += entity_ab as f64 * entity_mass;
    }

    Ok(mass)
}*/