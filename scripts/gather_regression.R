source(file.path('scripts', 'setup_regression.R'))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) == 0) {
  resultDirRaw = 'results/exome/20180716_082446'
} else {
  resultDirRaw = cmdArgs[1]}

############################################################

resultDirBase = sub(pattern = '^([^_]*_[^_]*)_.*$', replacement = '\\1',
                    x = basename(resultDirRaw))

resultDir = file.path(dirname(resultDirRaw), resultDirBase)
subDir = 'task_outputs'
dir.create(file.path(resultDir, subDir), recursive = TRUE)

taskDirs = intersect(list.files(dirname(resultDirRaw),
                                pattern = paste0(resultDirBase, '_.*'),
                                full.names = TRUE, include.dirs = TRUE),
                     list.dirs(dirname(resultDirRaw)))
taskFiles = c('progress_cox.txt', 'progress_logistic.txt', 'workspace.Rdata')

gwasMetadataList = gatherTaskResults(taskDirs, taskFiles, resultDir, subDir)

############################################################

incompleteIdx = sapply(gwasMetadataList, is.character)

gwasMetadata = rbindlist(gwasMetadataList[!incompleteIdx])
setkeyv(gwasMetadata, 'phecode')
write_tsv(gwasMetadata, file.path(resultDir, 'gwas_metadata.tsv'))

write_tsv(data.table(taskDirs[incompleteIdx]),
          file.path(resultDir, 'incomplete_tasks.tsv'), col_names = FALSE)

unlink(taskDirs[!incompleteIdx], recursive = TRUE)
