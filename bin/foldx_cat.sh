#!/usr/bin/env bash
# Script to combine FoldX output from multiple runs

PDB_ID=$1

# Process Average
head -n 1 Average_0_$PDB_ID.fxout > Average_$PDB_ID.fxout
tail -n +2 -q Average_*_$PDB_ID.fxout >> Average_$PDB_ID.fxout
#rm Average_*_$PDB_ID.fxout

# Process Dif
head -n 1 Dif_0_$PDB_ID.fxout > Dif_$PDB_ID.fxout
tail -n +2 -q Dif_*_$PDB_ID.fxout >> Dif_$PDB_ID.fxout
#rm Dif_*_$PDB_ID.fxout

# Process Raw
head -n 1 Raw_0_$PDB_ID.fxout > Raw_$PDB_ID.fxout
tail -n +2 -q Raw_*_$PDB_ID.fxout >> Raw_$PDB_ID.fxout
#rm Raw_*_$PDB_ID.fxout

# Process PdbList
head -n 1 PdbList_0_$PDB_ID.fxout > PdbList_$PDB_ID.fxout
tail -n +2 -q PdbList_*_$PDB_ID.fxout >> PdbList_$PDB_ID.fxout
#rm PdbList_*_$PDB_ID.fxout
