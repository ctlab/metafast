@echo off
SETLOCAL ENABLEDELAYEDEXPANSION

set this=%0
set mem=
set jvm_opts=
set params=

for %%x in (%*) do (
    set y=%%x

    if "!y!" == "-m" (
        set nextMem=1
    ) else if "!y!" == "--memory" (
        set nextMem=1
    ) else if "!nextMem!" == "1" (
        set mem=-Xmx!y! -Xms!y!
        set nextMem=0
    ) else if "!y!" == "-ea" (
        set jvm_opts=!jvm_opts! -ea
    ) else if "!y!" == "--enable-assertions" (
        set jvm_opts=!jvm_opts! -ea
    ) else if "!y:~0,2!" == "-X" (
        set jvm_opts=!jvm_opts! !y!
    ) else if "!y:~0,10!" == "-agentlib:" (
        set jvm_opts=!jvm_opts! !y!
    ) else (
        set params=!params! !y!
    )
)

if "!mem!" == "" (
    for /f "skip=1" %%p in ('wmic os get FreePhysicalMemory') do (
        if "!mem!" == "" set mem=%%p
    )
    set /a mem=mem/1024*90/100

    set JVM_64bit=0
    java -version 2>java_version
    for /f "skip=2 tokens=3" %%i in (java_version) do (
        if "%%i" == "64-Bit" (
            set JVM_64bit=1
        )
    )
    erase java_version
    if "!JVM_64bit!" == "0" (
        if !mem! gtr 1400 (
            set mem=1400
        )
    )

    set mem=-Xmx!mem!M -Xms!mem!M
)

java -Duser.language=en -Duser.country=US -Xss24M -XX:NewRatio=9 !mem! !jvm_opts! -jar !this! !params!
exit /b


