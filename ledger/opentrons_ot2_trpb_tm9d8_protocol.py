# Opentrons OT-2 Automated Mutagenesis Assembly Protocol
# Generated automatically by VeritasBio Wet-Lab Translator on 2026-10-28T01:09:49.489732+00:00
# Target: Tryptophan synthase beta chain, engineered parent Tm9D8* Variant Library (96-well format)

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
    "TRYP_S228A:A233I",
    "TRYP_A233I",
    "TRYP_G229Q",
    "TRYP_G230E",
    "TRYP_V227Y",
    "TRYP_G230I",
    "TRYP_N232C",
    "TRYP_G229F",
    "TRYP_F184C",
    "TRYP_F184I"
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
