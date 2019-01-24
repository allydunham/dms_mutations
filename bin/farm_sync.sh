#!/usr/bin/env bash
# Script to sync data and meta folders between farm and local folder
# Currently does not delete things - must be done manually for safety

#Colours for printf
green=$(tput setaf 2)
normal=$(tput sgr0)

## Sync To Farm
# Data folder
printf "%s\n%s\n%s\n" "${green}Rsyncing mutations project${normal}" "${green}Local -> Farm${normal}" "${green}File List: Data/${normal}"
rsync -vauh --exclude-from $HOME/.rsync_exclude  --dry-run $HOME/Projects/mutations/data/ ebi:/nfs/research1/beltrao/ally/mutations/data

read -p "Transfer? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   rsync -auh --info=progress2 --exclude-from $HOME/.rsync_exclude $HOME/Projects/mutations/data/ ebi:/nfs/research1/beltrao/ally/mutations/data
fi

# Meta folder
printf "\n%s\n" "${green}File List: Meta/${normal}"
rsync -vauh --exclude-from $HOME/.rsync_exclude --exclude=*.md --dry-run $HOME/Projects/mutations/meta/ ebi:/nfs/research1/beltrao/ally/mutations/meta

read -p "Transfer? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   rsync -auh --info=progress2 --exclude-from $HOME/.rsync_exclude --exclude=*.md $HOME/Projects/mutations/meta/ ebi:/nfs/research1/beltrao/ally/mutations/meta
fi

# Figures folder
printf "\n%s\n" "${green}File List: Figures/${normal}"
rsync -vauh --exclude-from $HOME/.rsync_exclude --exclude=*.md --dry-run $HOME/Projects/mutations/figures/ ebi:/nfs/research1/beltrao/ally/mutations/figures

read -p "Transfer? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   rsync -auh --info=progress2 --exclude-from $HOME/.rsync_exclude --exclude=*.md $HOME/Projects/mutations/figures/ ebi:/nfs/research1/beltrao/ally/mutations/figures
fi


## Sync From Farm
printf "\n%s\n" "${green}Farm -> Local${normal}"
# Data Folder
printf "%s\n" "${green}File List: Data${normal}"
rsync -vauh --exclude-from $HOME/.rsync_exclude  --dry-run ebi:/nfs/research1/beltrao/ally/mutations/data/ $HOME/Projects/mutations/data

read -p "Transfer? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   rsync -auh --info=progress2 --exclude-from $HOME/.rsync_exclude ebi:/nfs/research1/beltrao/ally/mutations/data/ $HOME/Projects/mutations/data
fi

# Meta Folder
printf "\n\n%s\n" "${green}File List: Meta${normal}"
rsync -vauh --exclude-from $HOME/.rsync_exclude --exclude=*.md --dry-run ebi:/nfs/research1/beltrao/ally/mutations/meta/ $HOME/Projects/mutations/meta

read -p "Transfer? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   rsync -auh --info=progress2 --exclude=*.md --exclude-from $HOME/.rsync_exclude ebi:/nfs/research1/beltrao/ally/mutations/meta/ $HOME/Projects/mutations/meta
fi

# Figures Folder
printf "\n\n%s\n" "${green}File List: Figures${normal}"
rsync -vauh --exclude-from $HOME/.rsync_exclude --exclude=*.md --dry-run ebi:/nfs/research1/beltrao/ally/mutations/figures/ $HOME/Projects/mutations/figures

read -p "Transfer? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   rsync -auh --info=progress2 --exclude=*.md --exclude-from $HOME/.rsync_exclude ebi:/nfs/research1/beltrao/ally/mutations/figures/ $HOME/Projects/mutations/figures
fi