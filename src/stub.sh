#!/usr/bin/env bash

jvm_opts=()
params=()

while [[ $# -ne 0 ]]
do
    if [[ "$1" == -m || "$1" == --memory ]]; then 
        mem=("${mem[@]}" "-Xmx$2" "-Xms$2")
        shift
    elif [[ "$1" == -ea || "$1" == --enable-assertions ]]; then 
        jvm_opts=("${jvm_opts[@]}" "-ea")
    elif [[ "$1" == -X* || "$1" == -agentlib:* ]]; then
        jvm_opts=("${jvm_opts[@]}" "$1")
    else
        params=("${params[@]}" "$1")
    fi
    shift
done

if [ ${#mem[@]} == 0 ]; then
    if which free > /dev/null 2>&1; then
        let mem=`free -m | grep "buffers/cache" | sed -r "s/\s+/ /g" | cut -d " " -f 4`
        mem=$(($mem * 95 / 100))M
        mem=("-Xmx$mem" "-Xms$mem")
    elif which vm_stat > /dev/null 2>&1; then
        inactive_mem=`vm_stat | grep -i inactive | grep -o "[[:digit:]]\+"`
        free_mem=`vm_stat | grep -i free | grep -o "[[:digit:]]\+"`
        mem=$((($inactive_mem + $free_mem) * 4096 / 1024 / 1024 * 95 / 100))M
        mem=("-Xmx$mem" "-Xms$mem")
    else 
        echo "WARN: Can't detect free memory, using java default"
    fi
fi

java -Xss24M -XX:NewRatio=9 "${mem[@]}" "${jvm_opts[@]}" -jar "$0" "${params[@]}"
exit


