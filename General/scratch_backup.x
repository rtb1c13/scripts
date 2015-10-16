#!/bin/bash

# Brief script to rsync whole scratch directory and keep backups in date-defined folders
# Source unknown, but not mine!

export date=`date "+%Y-%m-%dT%H:%M:%S"`
rsync -azP --link-dest=/media/Scratch2_backup/current /local/scratch/rtb1c13/ /media/Scratch2_backup/back-$date
rm -f /media/Scratch2_backup/current && ln -s /media/Scratch2_backup/back-$date /media/Scratch2_backup/current

