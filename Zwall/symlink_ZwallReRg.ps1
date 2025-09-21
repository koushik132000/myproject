#!/bin/bash
# fix_symlinks.sh

BASE="/home/laxmi/Desktop/koushik/myproject"
SHARED="../Zwall/ReCOMRgPBC"

# Define files as an array (space-separated, not comma-separated)
FILES=("3_EndtoEndZW.py" "4_COMZW.py" "5_RGZW.py" "10_Re_Zwall.py")

# List of project folders
PROJECTS=(
    "$BASE/LMloop_PBC_Zwall_seed1_Ea"
    "$BASE/LMloop_PBC_Zwall_seed2_Ea"
    "$BASE/LMloop_PBC_Zwall_seed1"
    "$BASE/LMloop_PBC_Zwall_seed2"
)

# Loop through projects and files
for project in "${PROJECTS[@]}"; do
  for file in "${FILES[@]}"; do
    link="$project/$file"
    target="$SHARED/$file"

    # Remove old file/link if exists
    [ -e "$link" ] && rm -f "$link"

    # Create symlink
    (cd "$project" && ln -s "$target" "$file")

    echo "Created symlink: $link -> $target"
  done
done
