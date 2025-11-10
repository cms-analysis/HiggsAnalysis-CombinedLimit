#!/bin/bash
# Version consistency checker for Combine releases
# This script verifies that version strings are consistent across the codebase

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

VERSION=$1

# Validate version format
if ! [[ $VERSION =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo -e "${RED}Error: Invalid version format${NC}"
    echo "Version must be in format vX.Y.Z (e.g., v10.3.2)"
    exit 1
fi

echo "Checking version consistency for ${VERSION}..."
echo ""

ERRORS=0

# Get script directory and repo root
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$( cd "$SCRIPT_DIR/.." && pwd )"

# Check 1: bin/combine.cpp
echo -n "Checking bin/combine.cpp... "
COMBINE_CPP="$REPO_ROOT/bin/combine.cpp"
if [ ! -f "$COMBINE_CPP" ]; then
    echo -e "${RED}FAIL${NC}"
    echo "  File not found: $COMBINE_CPP"
    ERRORS=$((ERRORS + 1))
else
    if grep -q "std::string combineTagString = \"${VERSION}\"" "$COMBINE_CPP"; then
        echo -e "${GREEN}OK${NC}"
    else
        echo -e "${RED}FAIL${NC}"
        CURRENT=$(grep "std::string combineTagString" "$COMBINE_CPP" | sed -E 's/.*"(.*)".*/\1/')
        echo "  Expected: std::string combineTagString = \"${VERSION}\""
        echo "  Found:    std::string combineTagString = \"${CURRENT}\""
        ERRORS=$((ERRORS + 1))
    fi
fi

# Check 2: docs/index.md - recommended tag line
echo -n "Checking docs/index.md (recommended tag)... "
DOCS_INDEX="$REPO_ROOT/docs/index.md"
if [ ! -f "$DOCS_INDEX" ]; then
    echo -e "${RED}FAIL${NC}"
    echo "  File not found: $DOCS_INDEX"
    ERRORS=$((ERRORS + 1))
else
    # Check for the main v10 recommended tag
    if grep -q "Currently, the recommended tag is \*\*${VERSION}\*\*" "$DOCS_INDEX"; then
        echo -e "${GREEN}OK${NC}"
    else
        echo -e "${RED}FAIL${NC}"
        CURRENT=$(grep "Currently, the recommended tag is" "$DOCS_INDEX" | head -1 | sed -E 's/.*\*\*([v0-9.]+)\*\*.*/\1/')
        echo "  Expected: Currently, the recommended tag is **${VERSION}**"
        echo "  Found:    Currently, the recommended tag is **${CURRENT}**"
        ERRORS=$((ERRORS + 1))
    fi
fi

# Check 3: docs/index.md - git clone branch
echo -n "Checking docs/index.md (git clone branch)... "
if grep -q "\-\-branch ${VERSION}" "$DOCS_INDEX"; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${RED}FAIL${NC}"
    CURRENT=$(grep "\-\-branch v" "$DOCS_INDEX" | head -1 | sed -E 's/.*--branch ([v0-9.]+).*/\1/')
    echo "  Expected: --branch ${VERSION}"
    echo "  Found:    --branch ${CURRENT}"
    ERRORS=$((ERRORS + 1))
fi

# Check 4: docs/index.md - release notes link
echo -n "Checking docs/index.md (release notes link)... "
if grep -q "releases/tag/${VERSION}" "$DOCS_INDEX"; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${YELLOW}WARNING${NC}"
    echo "  Release notes link not found: releases/tag/${VERSION}"
    echo "  This is expected if the release hasn't been created yet"
fi

echo ""
if [ $ERRORS -eq 0 ]; then
    echo -e "${GREEN}✓ All version checks passed!${NC}"
    echo ""
    echo "Next steps:"
    echo "  1. Commit the version updates:"
    echo "     git add bin/combine.cpp docs/index.md"
    echo "     git commit -m \"Update version to ${VERSION}\""
    echo "  2. Create and push the tag:"
    echo "     git tag -a ${VERSION} -m \"Release ${VERSION}\""
    echo "     git push origin main"
    echo "     git push origin ${VERSION}"
    exit 0
else
    echo -e "${RED}✗ Found ${ERRORS} error(s)${NC}"
    echo ""
    echo "Please fix the version inconsistencies and run this script again."
    exit 1
fi
