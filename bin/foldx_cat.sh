#!/usr/bin/env bash
# Script to combine FoldX output from multiple runs

PDB_ID=$1

# Process Average
head -n 9 Average_1_${PDB_ID}_repair.fxout > Average_${PDB_ID}_repair.fxout
tail -n +10 -q Average_*_${PDB_ID}_repair.fxout >> Average_$PDB_ID.fxout
#rm Average_*_${PDB_ID}_repair.fxout

# Process Dif
head -n 9 Dif_1_${PDB_ID}_repair.fxout > Dif_${PDB_ID}_repair.fxout
tail -n +10 -q Dif_*_${PDB_ID}_repair.fxout >> Dif_${PDB_ID}_repair.fxout
#rm Dif_*_${PDB_ID}_repair.fxout

# Process Raw
head -n 9 Raw_1_${PDB_ID}_repair.fxout > Raw_${PDB_ID}_repair.fxout
tail -n +10 -q Raw_*_${PDB_ID}_repair.fxout >> Raw_${PDB_ID}_repair.fxout
#rm Raw_*_${PDB_ID}_repair.fxout

# Process PdbList
cat PdbList_*_${PDB_ID}_repair.fxout >> PdbList_${PDB_ID}_repair.fxout
#rm PdbList_*_${PDB_ID}_repair.fxout
