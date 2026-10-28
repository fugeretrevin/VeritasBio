# Opentrons OT-2 Automated Mutagenesis Assembly Protocol
# Generated automatically by VeritasBio Wet-Lab Translator on 2026-10-27T17:16:49.013569+00:00
# Target: Aequorea victoria GFP Variant Library (96-well format)

from opentrons import protocol_api

metadata = {
    "protocolName": "VeritasBio High-Throughput Mutagenesis Assembly",
    "author": "VeritasBio Robotics Engine",
    "description": "Automated distribution of Gibson master mix, mutant primers, and linear template.",
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
    "GFP_F64L:Y66A:T203L",
    "GFP_F64L:Y66A:A206L",
    "GFP_F64L:Y66T:A206L",
    "GFP_F64L:S65A:Y66A",
    "GFP_F64L:Y66A",
    "GFP_F64A:S65K:T203L",
    "GFP_S65K:V163Y:T203L",
    "GFP_S65K:Y145L:T203L",
    "GFP_Y66T:S147T:A206L",
    "GFP_Y66A:A206L"
]

    protocol.comment(f"Starting VeritasBio automated assembly for {len(variants)} variants...")

    # Step 1: Distribute 10 uL 2X Gibson Assembly Master Mix to all active wells
    master_mix_well = reagent_reservoir.wells()[0]
    p300.pick_up_tip()
    for i in range(0, len(variants), 8):
        col_index = i // 8
        if col_index < 12:
            dest_col = destination_plate.columns()[col_index]
            p300.transfer(10, master_mix_well, dest_col, new_tip="never")
    p300.drop_tip()

    # Step 2: Transfer 2 uL of unique mutated primer / cDNA fragment to each destination well
    for idx, var_name in enumerate(variants):
        if idx >= 96:
            break
        src_well = primer_plate.wells()[idx]
        dst_well = destination_plate.wells()[idx]
        
        p20.pick_up_tip()
        p20.transfer(2.0, src_well, dst_well, mix_after=(3, 5), new_tip="never")
        p20.drop_tip()

    protocol.comment("Assembly complete. Transfer plate to thermocycler: 50 deg C for 60 minutes.")
