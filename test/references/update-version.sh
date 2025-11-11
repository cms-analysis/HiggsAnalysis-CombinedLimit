#!/bin/bash
# Update version strings in test reference files
# This script updates all *.out files in test/references/ with the new version

set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if version argument is provided
if [ $# -eq 0 ]; then
    echo -e "${RED}Error: No version specified${NC}"
    echo "Usage: $0 vX.Y.Z"
    echo "Example: $0 v10.3.2"
    exit 1
fi

NEW_VERSION=$1

# Validate version format
if ! [[ $NEW_VERSION =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo -e "${RED}Error: Invalid version format${NC}"
    echo "Version must be in format vX.Y.Z (e.g., v10.3.2)"
    exit 1
fi

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Updating test reference files to version ${NEW_VERSION}..."
echo ""

# Find current version in files
CURRENT_VERSION=$(grep -h "<<< v" "$SCRIPT_DIR"/*.out | head -1 | sed -E 's/.*<<< (v[0-9]+\.[0-9]+\.[0-9]+) >>>.*/\1/')

if [ -z "$CURRENT_VERSION" ]; then
    echo -e "${RED}Error: Could not find current version in test files${NC}"
    exit 1
fi

echo "Current version: ${CURRENT_VERSION}"
echo "New version: ${NEW_VERSION}"
echo ""

if [ "$CURRENT_VERSION" = "$NEW_VERSION" ]; then
    echo -e "${YELLOW}Warning: Version is already ${NEW_VERSION}${NC}"
    echo "No changes needed."
    exit 0
fi

# Update all .out files
UPDATED_COUNT=0
for file in "$SCRIPT_DIR"/*.out; do
    if [ -f "$file" ]; then
        if grep -q "<<< $CURRENT_VERSION >>>" "$file"; then
            sed -i.bak "s/<<< $CURRENT_VERSION >>>/<<< $NEW_VERSION >>>/g" "$file"
            rm "${file}.bak"
            echo -e "${GREEN}✓${NC} Updated: $(basename "$file")"
            UPDATED_COUNT=$((UPDATED_COUNT + 1))
        fi
    fi
done

echo ""
if [ $UPDATED_COUNT -eq 0 ]; then
    echo -e "${YELLOW}No files updated${NC}"
else
    echo -e "${GREEN}Successfully updated $UPDATED_COUNT file(s)${NC}"
fi
