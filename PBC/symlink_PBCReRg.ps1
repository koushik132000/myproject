# ============================
# setup_links.ps1
# Creates RELATIVE symbolic links for shared files in multiple projects
# Run this script in PowerShell (Developer Mode ON or as Administrator)
# ============================

# Relative path from project folder to shared folder
$SharedFolder = "..\ReCOMRgPBC"

# List of projects where the symlinks should be created
$Projects = @(
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\PBC\LMloopPBC_ideal_S1",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\PBC\LMloopPBC_ideal_S2"
)

# List of shared files
$Files = @(
    "3_EndtoEndPBC.py"
)

# Loop through projects and create symlinks
foreach ($project in $Projects) {
    foreach ($file in $Files) {
        $target = Join-Path $SharedFolder $file   # relative target
        $link   = Join-Path $project $file

        # Remove old file/link if it exists
        if (Test-Path $link) {
            Remove-Item $link -Force
        }

        # Create symbolic link inside project folder
        Push-Location $project
        New-Item -ItemType SymbolicLink -Path $file -Target $target | Out-Null
        Pop-Location

        Write-Host "Created symlink: $link -> $target"
    }
}
