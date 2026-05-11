$ErrorActionPreference = "Stop"

Write-Host "Initializing git repository..."
git init
git branch -M main
git remote add origin https://github.com/SubhadipJana1409/MultiOmicsBridge.git

Write-Host "Fetching list of files..."
# Get all files in the directory and subdirectories, excluding the .git folder
$files = Get-ChildItem -File -Recurse | Where-Object { $_.FullName -notmatch '\\\.git\\' }

$count = 0
foreach ($file in $files) {
    # Get relative path for git commands and cleaner commit messages
    $relativePath = Resolve-Path -Relative -Path $file.FullName
    $cleanPath = $relativePath -replace '^\.\\', ''
    $cleanPath = $cleanPath -replace '\\', '/'
    
    Write-Host "Adding and committing: $cleanPath"
    
    git add "`"$cleanPath`""
    # Only commit if there are changes (in case the file was already committed)
    $status = git status --porcelain "`"$cleanPath`""
    if ($status) {
        git commit -m "Add $cleanPath"
        $count++
    }
}

Write-Host "Successfully created $count commits."
Write-Host "You can now push the changes by running: git push -u origin main"
