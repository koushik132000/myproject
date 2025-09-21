# ============================
# setup_links.ps1
# Creates symbolic links for shared files in multiple projects
# Run this script in PowerShell (as Administrator)
# ============================

# Path to the folder containing the REAL shared files
$SharedPath = "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Zwall\ReCOMRgZwall"

# List of projects where the symlinks should be created
$Projects = @(
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Zwall\LMloop_PBC_Zwall_seed2")

# List of shared files
$Files = @(
    "6_LMloopPBC_Zwall.c"
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
