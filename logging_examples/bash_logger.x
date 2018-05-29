#!/bin/bash

# A (very) simple bash logger
# Run as 'bash_logger.x 2> logger_errfile.txt' to redirect stderr to file and catch error messages

DATETIME="date"
LOGFILE="logger_logfile.txt"
INFOFILE="logger_infofile.txt"

do_log() { echo `$DATETIME`" : $1" >> $LOGFILE; }
do_info() { echo `$DATETIME`" : $1" >> $INFOFILE; }

# Script starts here
do_info "Command line = $0 $@"
do_info "Script starting"
do_log "Command 1 running"
sleep 5
do_log "Command 1 success, errcode $?"
# Second command will fail
do_log "Command 2 running"
. nosuchfile.txt
do_log "Command 2 success, errcode $?"

