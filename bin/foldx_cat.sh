#!/usr/bin/env bash
# Script to combine FoldX output from multiple runs

PDB_ID=$1

# Process Average
head -n 9 Average_0_$PDB_ID.fxout > Average_$PDB_ID.fxout
tail -n +10 -q Average_*_$PDB_ID.fxout >> Average_$PDB_ID.fxout
#rm Average_*_$PDB_ID.fxout

# Process Dif
head -n 9 Dif_0_$PDB_ID.fxout > Dif_$PDB_ID.fxout
tail -n +10 -q Dif_*_$PDB_ID.fxout >> Dif_$PDB_ID.fxout
#rm Dif_*_$PDB_ID.fxout

# Process Raw
head -n 9 Raw_0_$PDB_ID.fxout > Raw_$PDB_ID.fxout
tail -n +10 -q Raw_*_$PDB_ID.fxout >> Raw_$PDB_ID.fxout
#rm Raw_*_$PDB_ID.fxout

# Process PdbList
cat PdbList_*_$PDB_ID.fxout >> PdbList_$PDB_ID.fxout
#rm PdbList_*_$PDB_ID.fxout
