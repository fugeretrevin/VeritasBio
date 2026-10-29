# Opentrons OT-2 Automated Mutagenesis Assembly Protocol
# Generated automatically by PanoptesBio Wet-Lab Translator on 2026-10-28T05:24:30.638403+00:00
# Target: Tryptophan synthase beta chain, engineered parent Tm9D8* Variant Library (96-well format)
# Anti-Edge Evaporation Barrier: True

from opentrons import protocol_api

metadata = {
    "protocolName": "PanoptesBio High-Throughput Mutagenesis Assembly",
    "author": "PanoptesBio Robotics Engine",
    "description": "Automated distribution of evaporation barriers, Gibson master mix, and mutant primers.",
    "apiLevel": "2.14"
}

def run(protocol: protocol_api.ProtocolContext):
    # Labware setup
    tips_20ul = protocol.load_labware("opentrons_96_tiprack_20ul", 1)
    tips_300ul = protocol.load_labware("opentrons_96_tiprack_300ul", 2)
    destination_plate = protocol.load_labware("corning_96_wellplate_360ul_flat", 3, "Gibson Assembly Plate")
    reagent_reservoir = protocol.load_labware("nest_12_reservoir_15ml", 4, "Master Mix & Buffers")
    primer_plate = protocol.load_labware("biorad_96_wellplate_200ul_pcr", 5, "Mutagenesis Primer Plate")

    # Pipette setup
    p20 = protocol.load_instrument("p20_single_gen2", "right", tip_racks=[tips_20ul])
    p300 = protocol.load_instrument("p300_multi_gen2", "left", tip_racks=[tips_300ul])

    variants = [
    "TRYP_V183L:F184I:S228G",
    "TRYP_V183I:F184I:S228G",
    "TRYP_F184I:S228G",
    "TRYP_F184I:V227I:S228G",
    "TRYP_F184V:S228G",
    "TRYP_F184L:S228G",
    "TRYP_F184M:S228G",
    "TRYP_V183M:F184I:S228G",
    "TRYP_V183A:F184I:S228G",
    "TRYP_F184I:V227T:S228G"
]
    perimeter_barrier_wells = [
    "A1",
    "A10",
    "A11",
    "A12",
    "A2",
    "A3",
    "A4",
    "A5",
    "A6",
    "A7",
    "A8",
    "A9",
    "B1",
    "B12",
    "C1",
    "C12",
    "D1",
    "D12",
    "E1",
    "E12",
    "F1",
    "F12",
    "G1",
    "G12",
    "H1",
    "H10",
    "H11",
    "H12",
    "H2",
    "H3",
    "H4",
    "H5",
    "H6",
    "H7",
    "H8",
    "H9"
]

    protocol.comment(f"Starting PanoptesBio automated assembly for {len(variants)} variants...")

    # Step 1: Dispense 200 uL sterile water/buffer evaporation barriers to eliminate edge effects
    barrier_reservoir_well = reagent_reservoir.wells()[1]
    p20_single = p20
    if True:
        protocol.comment("Dispensing 200 uL evaporation barrier to 36 perimeter wells...")
        p300.pick_up_tip()
        for well_name in perimeter_barrier_wells:
            dest_w = destination_plate[well_name]
            p300.transfer(200.0, barrier_reservoir_well, dest_w, new_tip="never")
        p300.drop_tip()

    # Step 2: Distribute 10 uL 2X Gibson Assembly Master Mix to active interior experimental wells
    master_mix_well = reagent_reservoir.wells()[0]
    p300.pick_up_tip()
    # Fill interior columns (2 through 11)
    for col_idx in range(1, 11):
        dest_col = destination_plate.columns()[col_idx]
        p300.transfer(10, master_mix_well, dest_col, new_tip="never")
    p300.drop_tip()

    # Step 3: Transfer 2 uL of unique mutated primer / cDNA fragment to each interior destination well
    # Experimental core starts at B6 (after B2:WT, B3:Neg, B4:Pos, B5:Blank)
    interior_rows = ["B", "C", "D", "E", "F", "G"]
    interior_wells = [f"{r}{c}" for r in interior_rows for c in range(2, 12)]
    candidate_wells = interior_wells[4:]  # Reserve first 4 for WT, Neg, Pos, Blank

    for idx, var_name in enumerate(variants):
        if idx >= len(candidate_wells):
            break
        src_well = primer_plate.wells()[idx]
        dst_well = destination_plate[candidate_wells[idx]]

        p20.pick_up_tip()
        p20.transfer(2.0, src_well, dst_well, mix_after=(3, 5), new_tip="never")
        p20.drop_tip()

    protocol.comment("Assembly complete. Transfer plate to thermocycler: 50 deg C for 60 minutes.")
