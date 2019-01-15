#!/usr/bin/env bash
# Script grep within a table and maintain the header
cat <(head -n 1 $2) <(grep $1 $2)