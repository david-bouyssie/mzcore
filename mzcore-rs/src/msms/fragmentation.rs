
use anyhow::*;
use std::collections::HashMap;
//use std::slice::Iter;
//use std::str::Chars;
//use itertools::Itertools;

use crate::chemistry::constants::WATER_MONO_MASS;
use crate::chemistry::table::AminoAcidTable;
use crate::ms::utils::*;
use crate::msms::model::*;

#[derive(Clone, PartialEq, Debug)]
pub struct SimpleFragmentationRule {
    pub fragment_ion_type: FragmentIonType,
    pub fragment_ion_charge: i8, // = 1
}

#[derive(Clone, Debug)]
pub struct SimpleFragmentationConfig {
    pub frag_rules: Vec<SimpleFragmentationRule>,
    pub frag_rule_by_ion_series: HashMap<FragmentIonSeries,SimpleFragmentationRule>
}

impl SimpleFragmentationConfig {

    /*pub fn create_with_defaults() -> Result<SimpleFragmentationConfig> {
        create(1 as i8, ActivationType::HCD, MsAnalyzer::FTMS)
    }*/

    /*pub fn create(
        fragment_ion_charge: i8,
        activation_type: ActivationType,
        msn_analyzer: MsAnalyzer
    ) -> Result<SimpleFragmentationConfig> {

        use ActivationType::*;
        use FragmentIonSeries::*;
        use MsAnalyzer::*;

        let ion_series_list = match activation_type  {
            CID => vec!(b,b_NH3,b_H2O,y,y_NH3,y_H2O),
            ECD => vec!(c,y,z_p1,z_p2),
            ETD => if msn_analyzer == FTMS {vec!(c,y,z,z_p1,z_p2)} else {vec!(c,y,z_p1,z_p2)},
            HCD => vec!(a,a_NH3,a_H2O,b,b_NH3,b_H2O,y,y_NH3,y_H2O,ya,yb),
            PSD => vec!(a,a_NH3,a_H2O,b,b_NH3,b_H2O,y),
        };

        /*logger.debug(
            s"Fragmentation config for activationType=$activationType and msnAnalyzer=$msnAnalyzer will contain this list of fragment series: " + ionSeriesList.mkString(",")
        )*/

        // Map fragmentation rules by series
        let mut frag_rule_by_ion_series:  HashMap<FragmentIonSeries,SimpleFragmentationRule> = HashMap::new();

        let simple_frag_rules = ion_series_list.iter().map( |ion_series| {

            let frag_rule = SimpleFragmentationRule {
                fragment_ion_type: FragmentIonType { ion_series: (*ion_series).to_owned(), neutral_loss: None },
                fragment_ion_charge: fragment_ion_charge
            };

            let frag_ion_series = frag_rule.fragment_ion_type.ion_series;
            if frag_rule_by_ion_series.contains_key(&frag_ion_series) {
                bail!("Duplicated SimpleFragmentationRule found, the same ionSeries {} has been already defined",frag_ion_series)
            }

            frag_rule_by_ion_series.insert(frag_ion_series, &*frag_rule);

            frag_rule
        }).collect();

        Ok(SimpleFragmentationConfig {
            frag_rules: simple_frag_rules,
            frag_rule_by_ion_series: frag_rule_by_ion_series,
        })
    }*/

    /*pub fn  containsIonSeries(&self, ionSeries: FragmentIonSeries.Value) -> bool {
        frag_rule_by_ion_series.contains(ionSeries)
    }
    pub fn  getFragmentationRule(&self, ionSeries: FragmentIonSeries.Value) -> Option[SimpleFragmentationRule] {
    fragRuleByIonSeries.get(ionSeries)
    }
    pub fn  getRequiredFragmentIonCharge(&self, ionSeries: FragmentIonSeries.Value): Option[Int] = {
    fragRuleByIonSeries.get(ionSeries).map(_.fragmentIonCharge)
    }*/
}

// --- Compute m/z value of fragment ion series (it contains all ion_types into one vector)  --- //
pub fn compute_frag_series_mz_values(pep_seq: &str, ion_type: FragmentIonSeries, charge: i8, aa_table: &AminoAcidTable) -> Result<Vec<f64>> {

    let seq_len = pep_seq.len();
    let mut frag_series_mz_values = Vec::with_capacity(seq_len);

    let frag_series_mass_shift = WATER_MONO_MASS + get_ion_mono_mass_shift(ion_type);

    let aa_by_code1 = &aa_table.aa_by_code1;

    let pep_seq_as_chars = pep_seq.chars();

    // --- This function contains ion series directions itself --- //
    use FragmentIonSeriesDirection::*;
    match get_ion_series_direction(ion_type) {
        Forward => {

            let mut seq_mass = frag_series_mass_shift;

            // Forward ions loop
            for aa_as_char in pep_seq_as_chars.take(seq_len - 1) {
                let aa_opt = aa_by_code1.get(&aa_as_char);
                let aa = aa_opt.context(format!(
                    "can't find amino acid '{}' in the provided table",aa_as_char
                ))?;

                seq_mass += aa.mono_mass;

                // Convert mass to m/z value
                let ion_mz = mass_to_mz(seq_mass, charge as i32);

                frag_series_mz_values.push(ion_mz);
            }

            // Forward ions loop
            /*for i in 1..seq_len { // starts at one because pep_seq slicing range will be between 0 and i-1
                let pep_subset_seq = pep_seq_as_chars[0..i].iter();
                //let pep_subset_seq = pep_seq_as_chars.take(i);
                frag_series_mz_values.push(_calc_frag_ion_mz(pep_subset_seq, aa_table, frag_series_mass_shift, charge)?);
                //println!("{} {:.4}", fion_aa, bmz);
            }*/

            Ok(frag_series_mz_values)
        }
        Reverse => {
            let mut seq_mass = frag_series_mass_shift;

            // Reverse ions loop
            for aa_as_char in pep_seq_as_chars.rev().take(seq_len - 1) {
                let aa_opt = aa_by_code1.get(&aa_as_char);
                let aa = aa_opt.context(format!(
                    "can't find amino acid '{}' in the provided table",aa_as_char
                ))?;

                seq_mass += aa.mono_mass;

                // Convert mass to m/z value
                let ion_mz = mass_to_mz(seq_mass, charge as i32);

                frag_series_mz_values.push(ion_mz);
            }

            // Reverse ions loop
            /*for i in (1..seq_len).rev() {
                // create slice (subset reference) of peptide sequence
                let pep_subset_seq = &pep_seq[i .. seq_len];
                frag_series_mz_values.push(_calc_frag_ion_mz(pep_subset_seq, aa_table, frag_series_mass_shift, charge)?);
            }*/

            Ok(frag_series_mz_values)
        }
        Unspecified => {
            bail!("unsupported ion type")
        }
    }

}

pub fn change_frag_series_charge_state(frag_series: &TheoreticalFragmentIons, new_charge: i8) -> TheoreticalFragmentIons {
    let mut new_mz_values = Vec::with_capacity(frag_series.mz_values.len());

    for mz_value in frag_series.mz_values.iter() {
        let mass = mz_to_mass(*mz_value, frag_series.charge as i32);
        new_mz_values.push( mass_to_mz(mass, new_charge as i32) );
    }

    TheoreticalFragmentIons {
        ion_type: frag_series.ion_type,
        charge: new_charge,
        mz_values: new_mz_values
    }
}

// --- Calculate mz value of a given peptide sequence depending on ion type and charge state --- //
/*fn _calc_frag_ion_mz(pep_subset_seq: Iter<char>, aa_table: &AminoAcidTable, frag_mass_shift: f64, charge: i8) -> Result<f64> {

    // Calculate mass of peptide sequence subset
    let pep_subset_mass = pep_subset_seq.calc_aa_seq_mass(pep_subset_seq, aa_table, true)?;

    // Calculate mass of fragment ion
    let ion_mass = pep_subset_mass + frag_mass_shift;

    // Convert mass to m/z value
    let ion_mz = mass_to_mz(ion_mass, charge as i32);
    //println!("{} {:.4}", rion_aa, ymz);

    Ok(ion_mz)
}*/

#[derive(Clone, PartialEq, Debug)]
pub struct TheoreticalFragmentIons {
    pub ion_type: FragmentIonSeries,
    pub charge: i8,
    pub mz_values: Vec<f64>
}

pub type FragmentationTable = Vec<TheoreticalFragmentIons>;

pub fn compute_frag_table_without_mods(pep_seq: &str, ion_types: &[FragmentIonSeries], frag_ion_charges: &Vec<i8>, aa_table: &AminoAcidTable) -> Result<FragmentationTable> {

    let mut frag_table: FragmentationTable = Vec::with_capacity(ion_types.len());

    // For each charge state, ion_types should be recalculated so that "for loop of charge" into the "for loop of ion_type" //
    // --- (e.g. b+1 for charge=1 and b+1 for charge=2) --- //
    for ion_type in ion_types {

        for charge in frag_ion_charges {

            let mz_values_res = compute_frag_series_mz_values(pep_seq, *ion_type, *charge, aa_table);
            let mz_values = mz_values_res?;

            // add to fragmentation table a new column containing different mz values for considered ion type and charge state
            frag_table.push(TheoreticalFragmentIons {
                ion_type: *ion_type,
                charge: *charge,
                mz_values: mz_values
            });
        }
    }

    // FragmentationTable
    Ok(frag_table)
}


pub fn compute_frag_table_with_mods(
    pep_seq: &str,
    frag_table_without_mods: &FragmentationTable,
    located_mass_increments: &Vec<(usize,f64)> // (aa_pos, mass_increment) // TODO: add custom Type
) -> Result<Vec<TheoreticalFragmentIons>> {

    let pep_seq_len = pep_seq.len();

    // Initial forward/reverse mass increments
    let mut forward_mass_increments = vec![0.0; pep_seq_len];
    let mut reverse_mass_increments = vec![0.0; pep_seq_len];

    // Inject mass increments in created vectors
    for &(aa_pos,mass_increment) in located_mass_increments {
        if aa_pos < 1 || aa_pos > pep_seq_len {
            bail!("invalid amino acid position ({}) for peptide {} of length {} (mod mass = {})", aa_pos, pep_seq, pep_seq_len, mass_increment);
        }
        forward_mass_increments[aa_pos - 1] += mass_increment;
        reverse_mass_increments[pep_seq_len - aa_pos] += mass_increment;
    }

    let mut cur_mass_inc = 0.0;
    for (idx,forward_mass_increment) in forward_mass_increments.clone().into_iter().enumerate() {
        cur_mass_inc += forward_mass_increment;
        forward_mass_increments[idx] = cur_mass_inc;
    }

    cur_mass_inc = 0.0;
    for (idx,reverse_mass_increment) in reverse_mass_increments.clone().into_iter().enumerate() {
        cur_mass_inc += reverse_mass_increment;
        reverse_mass_increments[idx] = cur_mass_inc;
    }

    let forward_mass_increments_ref = &forward_mass_increments;
    let reverse_mass_increments_ref = &reverse_mass_increments;

    //use FragmentIonSeries::*;
    //let default_frag_table= compute_frag_table_without_mods(pep_seq, ion_types, frag_ion_charges, aa_table)?;
    //println!("frag_table_prec_3: {:?}",frag_table_prec_3)

    let mut updated_frag_table = Vec::with_capacity(frag_table_without_mods.len());
    for frag_series in frag_table_without_mods.iter() {

        let mut updated_frag_mz_values = Vec::with_capacity(frag_series.mz_values.len());
        let mut cloned_frag_series = frag_series.clone();

        let mass_increments = if is_ion_forward(frag_series.ion_type).unwrap() {
            forward_mass_increments_ref
        } else {
            reverse_mass_increments_ref
        };

        let frag_series_charge = cloned_frag_series.charge as f64;
        for (frag_idx, frag_mz) in cloned_frag_series.mz_values.into_iter().enumerate() {
            let  mass_increment = mass_increments[frag_idx] / frag_series_charge;
            updated_frag_mz_values.push(frag_mz + mass_increment);
        }

        cloned_frag_series.mz_values = updated_frag_mz_values;
        updated_frag_table.push(cloned_frag_series);
    }

    Ok(updated_frag_table)
}

pub fn compute_frag_table_from_mod_string(
    pep_seq: &str,
    pep_mods_str_opt: Option<&str>,
    ion_types: &[FragmentIonSeries],
    frag_ion_charges: &Vec<i8>,
    aa_table: &AminoAcidTable
) -> Result<FragmentationTable> {

    let frag_table_without_mods = compute_frag_table_without_mods(
        pep_seq,
        ion_types,
        frag_ion_charges,
        aa_table
    )?;

    let frag_table = if pep_mods_str_opt.is_none() {
        frag_table_without_mods
    } else {
        let pep_mods_str = pep_mods_str_opt.unwrap();
        let pep_mods = pep_mods_str.split(",");

        let mut located_mass_incs = Vec::new();
        for pep_mod in pep_mods {
            let pep_mod_parts: Vec<_> = pep_mod.splitn(2, "@").collect();
            let mod_mass: f64 = pep_mod_parts.first().unwrap().parse().unwrap();
            let mut mod_pos: i64 = pep_mod_parts.last().unwrap().parse().unwrap();
            if mod_pos == 0 {
                mod_pos += 1;
            } else if mod_pos == -1 {
                mod_pos = pep_seq.len() as i64;
            }

            located_mass_incs.push((mod_pos as usize, mod_mass));
        }

        compute_frag_table_with_mods(
            pep_seq,
            &frag_table_without_mods,
            &located_mass_incs,
        )?
    };

    Ok(frag_table)
}

pub fn compute_frag_table(
    pep_seq: &str,
    located_mass_increments: &Vec<(usize,f64)>,
    ion_types: &[FragmentIonSeries],
    frag_ion_charges: &Vec<i8>,
    aa_table: &AminoAcidTable
) -> Result<FragmentationTable> {

    let frag_table_without_mods = compute_frag_table_without_mods(
        pep_seq,
        ion_types,
        frag_ion_charges,
        aa_table
    )?;

    let frag_table = compute_frag_table_with_mods(
        pep_seq,
        &frag_table_without_mods,
        &located_mass_increments,
    )?;

    Ok(frag_table)
}