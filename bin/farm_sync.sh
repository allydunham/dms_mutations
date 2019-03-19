#!/usr/bin/env bash
# Script to sync data and meta folders between farm and local folder
# Currently does not delete things - must be done manually for safety

#Colours for printf
green=$(tput setaf 2)
normal=$(tput sgr0)

# Local/Farm Roots
local=$HOME/Projects/mutations
farm=ebi:/nfs/research1/beltrao/ally/mutations

# Sync function
syncr () {
   rsync -vauh --exclude-from $HOME/.rsync_exclude --exclude-from $local/rsync_exclude --dry-run $1 $2

   read -p "Transfer? " -n 1 -r
   echo
   if [[ $REPLY =~ ^[Yy]$ ]]
   then
      rsync -auh --exclude-from $HOME/.rsync_exclude --exclude-from $local/rsync_exclude $1 $2
   fi
}

## Sync To Farm
printf "%s\n%s\n%s\n" "${green}Rsyncing mutations project${normal}" "${green}Local -> Farm${normal}" "${green}File List: Data/${normal}"
syncr $local/data/ $farm/data

printf "\n%s\n" "${green}File List: Meta/${normal}"
syncr $local/meta/ $farm/meta

printf "\n%s\n" "${green}File List: Figures/${normal}"
syncr $local/figures/ $farm/figures

## Sync From Farm
printf "\n%s\n%s\n" "${green}Farm -> Local${normal}" "${green}File List: Data${normal}"
syncr $farm/data/ $local/data

printf "\n\n%s\n" "${green}File List: Meta${normal}"
syncr $farm/meta/ $local/meta

printf "\n\n%s\n" "${green}File List: Figures${normal}"
syncr $farm/figures/ $local/figures