<#
.SYNOPSIS
Extract searchable text from PDFs registered for Manuscript 3.

.DESCRIPTION
The original PDF remains authoritative. This script creates a local text cache
for targeted search and reading. It does not perform OCR, rasterise pages, copy
PDFs into the repository, or create evidence notes automatically.

.EXAMPLE
.\literature_ingest.ps1 -List

.EXAMPLE
.\literature_ingest.ps1 -Key Declerck2025,Janke2022

.EXAMPLE
.\literature_ingest.ps1 -AllSources
#>

[CmdletBinding(DefaultParameterSetName = 'Keys')]
param(
    [Parameter(ParameterSetName = 'Keys')]
    [string[]]$Key,

    [Parameter(ParameterSetName = 'All')]
    [switch]$AllSources,

    [Parameter(ParameterSetName = 'List')]
    [switch]$List,

    [string]$Registry,
    [string]$Cache,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$workflowRoot = Split-Path -Parent $PSScriptRoot
if (-not $Registry) {
    $Registry = Join-Path $workflowRoot 'knowledge_base\00_source_registry\paper_registry.csv'
}
if (-not $Cache) {
    $Cache = Join-Path $workflowRoot 'knowledge_base\02_source_cache\extracted_text'
}
$manifest = Join-Path $workflowRoot 'knowledge_base\02_source_cache\extraction_manifest.csv'

if (-not (Test-Path -LiteralPath $Registry -PathType Leaf)) {
    throw "Source registry not found: $Registry"
}

$sources = @(Import-Csv -LiteralPath $Registry -Encoding UTF8)
if ($sources.Count -eq 0) {
    throw "Source registry is empty: $Registry"
}

if ($List) {
    $sources | ForEach-Object {
        $sourcePath = if ($_.pdf_path) { $_.pdf_path } else { 'SOURCE NOT REGISTERED' }
        '{0}: {1} | {2}' -f $_.citation_key, $_.year, $sourcePath
    }
    exit 0
}

if ($AllSources) {
    $selected = $sources
}
else {
    if (-not $Key -or $Key.Count -eq 0) {
        throw 'Select -Key, -AllSources or -List. No extraction was performed.'
    }
    $unknown = @($Key | Where-Object { $_ -notin $sources.citation_key })
    if ($unknown.Count -gt 0) {
        throw ('Citation key not found: ' + ($unknown -join ', '))
    }
    $selected = @($sources | Where-Object { $_.citation_key -in $Key })
}

$pdftotext = Get-Command pdftotext -ErrorAction SilentlyContinue
if (-not $pdftotext) {
    throw 'pdftotext was not found on PATH.'
}

New-Item -ItemType Directory -Force -Path $Cache | Out-Null
$processedAt = [DateTimeOffset]::UtcNow.ToString('o')
$records = foreach ($source in $selected) {
    $outputPath = Join-Path $Cache ($source.citation_key + '.txt')
    $status = ''
    $characters = 0
    $message = ''

    if (-not $source.pdf_path) {
        $status = 'source_not_registered'
    }
    elseif (-not (Test-Path -LiteralPath $source.pdf_path -PathType Leaf)) {
        $status = 'source_missing'
    }
    elseif ((Test-Path -LiteralPath $outputPath -PathType Leaf) -and -not $Force) {
        $characters = (Get-Content -LiteralPath $outputPath -Raw -Encoding UTF8).Length
        $status = 'existing_cache'
        $message = 'Use -Force to replace.'
    }
    else {
        & $pdftotext.Source -layout -enc UTF-8 -- $source.pdf_path $outputPath
        if ($LASTEXITCODE -ne 0) {
            if (Test-Path -LiteralPath $outputPath) {
                Remove-Item -LiteralPath $outputPath -Force
            }
            $status = 'pdftotext_failed'
        }
        else {
            $characters = (Get-Content -LiteralPath $outputPath -Raw -Encoding UTF8).Trim().Length
            if ($characters -ge 500) {
                $status = 'extracted'
            }
            else {
                $status = 'needs_ocr_or_review'
                $message = 'Very little text was extracted.'
            }
        }
    }

    [pscustomobject]@{
        citation_key = $source.citation_key
        pdf_path = $source.pdf_path
        text_path = $outputPath
        status = $status
        characters = $characters
        processed_at = $processedAt
        message = $message
    }
}

$records | Export-Csv -LiteralPath $manifest -NoTypeInformation -Encoding UTF8
$records | ForEach-Object {
    '{0}: {1} ({2} characters)' -f $_.citation_key, $_.status, $_.characters
}

$failureStates = @(
    'source_not_registered',
    'source_missing',
    'pdftotext_failed',
    'needs_ocr_or_review'
)
if (@($records | Where-Object { $_.status -in $failureStates }).Count -gt 0) {
    exit 1
}
exit 0
