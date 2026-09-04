#define MyAppName "Al Fatigue"
#define MyAppVersion "0.1.0"
#define MyAppPublisher "Al Fatigue Project"
#define MyAppExeName "AlFatigue.exe"

[Setup]
AppId={{7A5E2F33-DB7D-4C27-9F0E-0A1F47100001}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={autopf}\AlFatigue
DefaultGroupName={#MyAppName}
OutputDir=..\..\..\dist\installer
OutputBaseFilename=AlFatigue-Setup-{#MyAppVersion}
Compression=lzma
SolidCompression=yes
PrivilegesRequired=lowest
ArchitecturesInstallIn64BitMode=x64compatible
UninstallDisplayIcon={app}\{#MyAppExeName}

[Files]
Source: "..\..\..\dist\AlFatigue\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"
Name: "{autodesktop}\{#MyAppName}"; Filename: "{app}\{#MyAppExeName}"

[Registry]
Root: HKCU; Subkey: "Software\Classes\AlFatigue.Project"; ValueType: string; ValueName: ""; ValueData: "Al Fatigue Simulation Project"; Flags: uninsdeletekey
Root: HKCU; Subkey: "Software\Classes\AlFatigue.Project\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\{#MyAppExeName}"" ""%1"""
Root: HKCU; Subkey: "Software\Classes\.ftgsim"; ValueType: string; ValueName: ""; ValueData: "AlFatigue.Project"

[Run]
Filename: "{app}\{#MyAppExeName}"; Description: "Launch {#MyAppName}"; Flags: nowait postinstall skipifsilent
