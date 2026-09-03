param(
    [Parameter(Mandatory = $true)]
    [string]$Executable
)

$ErrorActionPreference = "Stop"
$exe = (Resolve-Path $Executable).Path
$progId = "AlFatigue.Project"
$root = "HKCU:\Software\Classes"

New-Item -Path "$root\$progId\shell\open\command" -Force | Out-Null
Set-ItemProperty -Path "$root\$progId" -Name '(Default)' -Value 'Al Fatigue Simulation Project'
Set-ItemProperty -Path "$root\$progId\shell\open\command" -Name '(Default)' -Value ('"{0}" "%1"' -f $exe)

New-Item -Path "$root\.ftgsim" -Force | Out-Null
Set-ItemProperty -Path "$root\.ftgsim" -Name '(Default)' -Value $progId

Write-Host ".ftgsim is now associated with $exe for the current user."
