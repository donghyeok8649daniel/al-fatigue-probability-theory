param([Parameter(Mandatory=$true)][string]$Pythonw)
$ErrorActionPreference = 'Stop'
$repo = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Pythonw = (Resolve-Path $Pythonw).Path
$install = Join-Path $env:LOCALAPPDATA 'AlFatigueProbabilityUI'
New-Item -ItemType Directory -Path $install -Force | Out-Null
$exe = Join-Path $install 'AlFatigueProbability.exe'
if (Test-Path $exe) { throw 'Launcher already exists; inspect it before replacing.' }
$compiler = Join-Path $env:WINDIR 'Microsoft.NET/Framework64/v4.0.30319/csc.exe'
& $compiler /nologo /target:winexe /reference:System.Windows.Forms.dll "/out:$exe" (Join-Path $PSScriptRoot 'windows_launcher.cs')
if ($LASTEXITCODE -ne 0) { throw 'Launcher compilation failed' }
[IO.File]::WriteAllLines((Join-Path $install 'launcher.paths'), @($Pythonw,$repo), [Text.UTF8Encoding]::new($false))
$desktop = [Environment]::GetFolderPath('Desktop')
if ([string]::IsNullOrWhiteSpace($desktop) -or !(Test-Path $desktop)) { throw 'Windows Desktop path unavailable' }
$link = Join-Path $desktop 'Al Fatigue UI.lnk'
if (Test-Path $link) { throw 'Shortcut already exists; inspect it before replacing.' }
$shell = New-Object -ComObject WScript.Shell
$shortcut = $shell.CreateShortcut($link)
$shortcut.TargetPath = $exe
$shortcut.WorkingDirectory = $repo
$shortcut.Description = 'Al Fatigue probability PDE - current research UI'
$shortcut.Save()
$prog = 'HKCU:/Software/Classes/AlFatigueProbability.Project'
New-Item -Path "$prog/shell/open/command" -Force | Out-Null
Set-Item -Path $prog -Value 'Al Fatigue probability project'
Set-Item -Path "$prog/shell/open/command" -Value ('"'+$exe+'" "%1"')
New-Item -Path 'HKCU:/Software/Classes/.ftgsim/OpenWithProgids' -Force | Out-Null
New-ItemProperty -Path 'HKCU:/Software/Classes/.ftgsim/OpenWithProgids' -Name 'AlFatigueProbability.Project' -Value '' -PropertyType String -Force | Out-Null
# Preserve an existing default and Windows UserChoice. Register Open With in all cases.
$ext = Get-Item 'HKCU:/Software/Classes/.ftgsim'
if ([string]::IsNullOrWhiteSpace($ext.GetValue(''))) { Set-Item $ext.PSPath -Value 'AlFatigueProbability.Project' }
Write-Output "Launcher: $exe"
Write-Output "Shortcut: $link"
Write-Output 'Registered .ftgsim Open With; existing default/UserChoice retained.'
