#!/usr/bin/env bash
# Script to sync folders between farm and local folder
# Currently does not delete things - must be done manually for safety
# Syncs folders passeed as arguments, or all folders if no args passed

# Colours for printf
green=$(tput setaf 2)
magenta=$(tput setaf 5)
bold=$(tput bold)
normal=$(tput sgr0)

# Local/Farm Roots
local=$HOME/Projects/mutations
farm=ebi:/nfs/research1/beltrao/ally/mutations

# Sync function
syncr () {
   rsync -vauh --exclude-from "$HOME/.rsync_exclude" --exclude-from "$local/rsync_exclude" --dry-run "$1" "$2"

   read -p "Transfer? " -n 1 -r
   echo
   if [[ $REPLY =~ ^[Yy]$ ]]
   then
      rsync -auh --exclude-from "$HOME/.rsync_exclude" --exclude-from "$local/rsync_exclude" "$1" "$2"
   fi
}

if [ $# -eq 0 ]
   then
   folders=( "data" "meta" "figures" )
else
   folders=( "$@" )
fi

printf "%s" "${magenta}${bold}Rsyncing mutations project${normal}"
for f in "${folders[@]}"
do
   printf "\n%s\n%s\n" "${green}${bold}Folder: $f${normal}" "${green}Local -> Farm${normal}"
   syncr "$local/$f/" "$farm/$f"
   printf "\n%s\n" "${green}Farm -> Local${normal}"
   syncr "$farm/$f/" "$local/$f"
done
