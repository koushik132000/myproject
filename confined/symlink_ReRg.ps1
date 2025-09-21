# ============================
# setup_links.ps1
# Creates symbolic links for shared files in multiple projects
# Run this script in PowerShell (as Administrator)
# ============================

# Path to the folder containing the REAL shared files
$SharedPath = "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\confined\ReCOMRgconfined"

# List of projects where the symlinks should be created
$Projects = @(
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\confined\LMloopMC_seed1",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\confined\LMloopMC_seed2")

# List of shared files
$Files = @(
    "3_EndtoEnd.py",
    "4_COM.py",
    "5_RG.py"
)

# Loop through projects and create symlinks
foreach ($project in $Projects) {
    foreach ($file in $Files) {
        $target = Join-Path $SharedPath $file
        $link   = Join-Path $project $file

        # Remove old file/link if it exists
        if (Test-Path $link) {
            Remove-Item $link -Force
        }

        # Create symbolic link
        New-Item -ItemType SymbolicLink -Path $link -Target $target | Out-Null
        Write-Host "Created symlink: $link -> $target"
    }
}
