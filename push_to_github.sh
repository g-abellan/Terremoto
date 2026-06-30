#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

echo "=================================================="
echo "GitHub Automation Script for Terremoto Project"
echo "=================================================="

# Check if git is installed
if ! command -v git &> /dev/null; then
    echo "Error: git is not installed on your system."
    exit 1
fi

# Initialize git repository if not already initialized
if [ ! -d ".git" ]; then
    echo "Initializing local Git repository..."
    git init
else
    echo "Git repository already initialized."
fi

# Add all files to staging
echo "Adding files to Git..."
git add .

# Check if there are changes to commit
if git diff --cached --quiet; then
    echo "No new changes to commit."
else
    echo "Committing files..."
    git commit -m "Initial commit: replication code, stabilized parameters, and manuscript draft"
fi

# Rename branch to main
echo "Setting main branch..."
git branch -M main

# Prompt for the GitHub URL
echo ""
echo "Please enter the HTTPS clone URL of your empty GitHub repository"
echo "Example: https://github.com/username/Terremoto.git"
read -p "Repository URL: " repo_url

if [ -z "$repo_url" ]; then
    echo "Error: Repository URL cannot be empty."
    exit 1
fi

# Remove existing origin if it exists
if git remote | grep -q "^origin$"; then
    echo "Updating existing remote origin..."
    git remote remove origin
fi

# Add remote origin
git remote add origin "$repo_url"

# Push to GitHub
echo "Pushing code to GitHub main branch..."
git push -u origin main

echo ""
echo "=================================================="
echo "Success! Your project has been uploaded to GitHub."
echo "=================================================="
