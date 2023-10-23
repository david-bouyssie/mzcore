///
/// Some parts of this file originates from [Sage](https://github.com/lazear/sage/blob/master/crates/sage/src/mass.rs)
/// Copyright (c) 2022 Michael Lazear
/// SPDX-License-Identifier: MIT
///
///
use anyhow::*;
use std::collections::HashMap;

pub fn parse_aa_composition(sequence: &str) -> Result<HashMap<char, f32>> {

    let seq_chars: Vec<char> = sequence.chars().filter(|c| !c.is_whitespace()).collect();
    let seq_len = seq_chars.len();

    // Count the AA occurrences
    let mut aa_count_by_char: HashMap<char, i32> = HashMap::new();

    let mut i = 0;
    while i < seq_len {
        let aa = seq_chars[i];
        let counter = aa_count_by_char.entry(aa).or_insert(0);
        *counter += 1;
        i += 1
    }

    let abundance_map: HashMap<char, f32> = aa_count_by_char.into_iter().map(|e| (e.0,e.1 as f32) ).collect();

    // Build the abundanceMap
    /*let mut abundance_map: HashMap<char, f32> = HashMap::new();

    for (aaCharAsInt,aaCount) in aa_count_by_char {

        let aaOpt = aaByCode1.get(aaCharAsInt)
        require(aaOpt.isDefined, s"amino acid ${aaCharAsInt.toChar} is missing in provided aaTable")

        abundance_map.insert(aaOpt.unwrap(), aaCount)
    }*/

    Ok(abundance_map)
}

/*pub fn parse_atom_composition(formula: &str, atom_table: AtomTable) -> Result<HashMap<char, f32>> {

    let atom_by_symbol = atom_table.atom_by_symbol;*/

pub fn parse_atom_composition(formula: &str) -> Result<HashMap<String, f32>> {

    let mut abundance_map: HashMap<String, f32> = HashMap::new();

    let formula_elements = formula.split(" ");

    for elem in formula_elements {

        let (atom_symbol, abundance) = if elem.contains('(') == false { (elem,1) }
        else {
            let elem_ab_parts: Vec<&str> = elem.split('(').collect();
            let elem_symbol = *elem_ab_parts.first().ok_or_else(|| anyhow!("no element symbol"))?;
            let elem_quant_str = *elem_ab_parts.last().ok_or_else(|| anyhow!("no element abundance"))?;
            let elem_quant: i32 = elem_quant_str.replace(')', "").parse()?;

            (elem_symbol, elem_quant)
        };

        /*if atom_symbol.len() != 1 {
            bail!("Invalid atom symbol '{}'",atom_symbol);
        }*/

        //let atom_opt = atom_by_symbol.get(atom_symbol);
        //let atom = atom_opt.ok_or_else(|| anyhow!("atom symbol {} is missing in provided atom_table",atom_symbol))?;
        //let atom = atom_symbol.chars().next().unwrap();

        abundance_map.insert(atom_symbol.to_string(), abundance as f32);
    }

    Ok(abundance_map)
}

pub fn sum_atom_compositions(atom_comp1: HashMap<String, f32>, atom_comp2: HashMap<String, f32>) -> HashMap<String, f32> {
    let mut new_atom_comp = atom_comp1.clone();
    for (atom,ab) in atom_comp2 {
        *new_atom_comp.entry(atom).or_insert(0.0) += ab;
    }

    new_atom_comp
}


/*
// --- Sage definitions --- //
use std::iter::Sum;

pub const fn composition(aa: u8) -> Composition {
    match aa {
        b'A' => Composition::new(3, 2, 0),
        b'R' => Composition::new(6, 2, 0),
        b'N' => Composition::new(4, 3, 0),
        b'D' => Composition::new(4, 4, 0),
        b'C' => Composition::new(3, 2, 1),
        b'E' => Composition::new(5, 4, 0),
        b'Q' => Composition::new(5, 3, 0),
        b'G' => Composition::new(2, 2, 0),
        b'H' => Composition::new(6, 2, 0),
        b'I' => Composition::new(6, 2, 0),
        b'L' => Composition::new(6, 2, 0),
        b'K' => Composition::new(6, 2, 0),
        b'M' => Composition::new(5, 2, 1),
        b'F' => Composition::new(9, 2, 0),
        b'P' => Composition::new(5, 2, 0),
        b'S' => Composition::new(3, 3, 0),
        b'T' => Composition::new(4, 3, 0),
        b'W' => Composition::new(11, 2, 0),
        b'Y' => Composition::new(9, 3, 0),
        b'V' => Composition::new(5, 2, 0),
        b'U' => Composition::new(3, 2, 0),
        b'O' => Composition::new(12, 3, 0),
        _ => Composition::new(0, 0, 0),
    }
}

pub struct Composition {
    pub carbon: u16,
    pub sulfur: u16,
}

impl Composition {
    pub const fn new(carbon: u16, _oxygen: u8, sulfur: u16) -> Self {
        Self { carbon, sulfur }
    }
}

impl Sum for Composition {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut comp = Composition::new(0, 0, 0);
        for i in iter {
            comp.carbon += i.carbon;
            // comp.oxygen += i.oxygen;
            comp.sulfur += i.sulfur;
        }
        comp
    }
}*/