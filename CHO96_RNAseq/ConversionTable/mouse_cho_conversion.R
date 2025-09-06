mcho <- fread('databases/ConversionTable/mouse_cho_biomart.tsv', na.strings = '')
mcho.parsed <- data.table()
mcho.parsed %<>% rbind(mcho[is.na(`Chinese hamster CriGri gene name`)&!is.na(`Chinese hamster CHOK1GS gene name`), .(mouseSymbol = `Gene name`, choSymbol = `Chinese hamster CHOK1GS gene name`)])

mcho.parsed %<>% rbind(mcho[!is.na(`Chinese hamster CriGri gene name`)&is.na(`Chinese hamster CHOK1GS gene name`), .(mouseSymbol = `Gene name`, choSymbol = `Chinese hamster CriGri gene name`)])

#for ambiguous mapping
mcho.parsed %<>% rbind(mcho[`Chinese hamster CriGri gene name`!=`Chinese hamster CHOK1GS gene name`] %>% 
                         melt(id.vars = 'Gene name') %>% 
                         .[, .(choSymbol = unique(value)),by = 'Gene name'] %>% setnames('Gene name', 'mouseSymbol')
)
mcho.parsed %<>% rbind(mcho[`Chinese hamster CriGri gene name`==`Chinese hamster CHOK1GS gene name`] %>% .[, .(mouseSymbol = `Gene name`, choSymbol = `Chinese hamster CHOK1GS gene name`)])
mcho.parsed %>% unique() %>% fwrite('databases/ConversionTable/mouse_cho_biomart_clean.tsv')
