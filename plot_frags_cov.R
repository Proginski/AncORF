# For each fragment obtain in the ancestrale reconstruction of de novo CDS, plot its maximum overlapp with another frag generated for the same ancestor
data <- read.table("/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SCER_SGD/max_covs.txt")
hist(data[,2])
