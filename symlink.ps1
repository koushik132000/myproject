# ============================
# setup_links.ps1
# Creates symbolic links for shared files in multiple projects
# Run this script in PowerShell (as Administrator)
# ============================

# Path to the folder containing the REAL shared files
$SharedPath = "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Rootmeanstderror"

# List of projects where the symlinks should be created
$Projects = @(
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Zwall\LMloop_PBC_Zwall_seed1_Ea",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Zwall\LMloop_PBC_Zwall_seed2_Ea",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Zwall\LMloop_PBC_Zwall_seed1",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\Zwall\LMloop_PBC_Zwall_seed2",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\PBC\LMloopPBC_ideal_S1",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\PBC\LMloopPBC_ideal_S2",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\PBC\LMloopPBC_MC_seed1",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\PBC\LMloopPBC_MC_seed2",
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\confined\LMloopMC_seed1"
    "C:\Users\Koushik Sai\OneDrive\Desktop\polymers\confined\LMloopMC_seed2")

# List of shared files
$Files = @(
    "8_RootMeanStdError.py",
    "9_RgRootMeanStdError.py",
    "10_Rgx_Rgy_Rgz.py"
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
