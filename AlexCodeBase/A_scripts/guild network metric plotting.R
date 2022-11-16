# plotting network stats through time
# Toarcian

stats <- read.csv("./data/metrics_time.csv")

guilds <- stats[grep(TRUE,stats[,"resolution"] == "guild"),]

########## STRUCTURE

# richness
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$taxa, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0,60), xlab="interval", ylab="Number of guilds")
  axis(2)
  axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# connectance
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$connectance, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.1,0.2), xlab="interval", ylab="Connectance (C)")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# max TL
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$max_tl, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(3,4), xlab="interval", ylab="Maximum trophic level (maxTL)")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# generality
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$generality, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.15,0.25), xlab="interval", ylab="Generality")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# vulnerability
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$vulnerability, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.1,0.15), xlab="interval", ylab="Vulnerability")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

########## multiplot
par(mfrow=c(3,2))
# richness
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$taxa, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0,60), xlab="interval", ylab="Number of guilds")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 60, "A")

# connectance
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$connectance, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.1,0.2), xlab="interval", ylab="Connectance (C)")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 0.2, "B")


# max TL
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$max_tl, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(3,4), xlab="interval", ylab="Maximum trophic level (maxTL)")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 4, "C")

# generality
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$generality, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.15,0.25), xlab="interval", ylab="Generality")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 0.25, "D")

# vulnerability
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$vulnerability, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.1,0.15), xlab="interval", ylab="Vulnerability")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 0.15, "E")

par(mfrow=c(1,1))

############ MOTIFS

# motif S1 (no. of linear chains/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s1, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.2,0.4), xlab="interval", ylab="S1: number of linear chains")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# motif S2 (no. of omnivory motifs/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s2, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.1,0.4), xlab="interval", ylab="S2: omnivory")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# motif S4 (no. of apparent competition motifs/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s4, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.35,1), xlab="interval", ylab="S4: apparent competition")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

# motif S5 (no. of direct competition motifs/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s5, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.15,0.35), xlab="interval", ylab="S5: direct competition")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()

########## multiplot
par(mfrow=c(2,2))

# motif S1 (no. of linear chains/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s1, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.2,0.4), xlab="interval", ylab="S1: number of linear chains")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 0.4, "A")

# motif S2 (no. of omnivory motifs/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s2, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.1,0.4), xlab="interval", ylab="S2: omnivory")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 0.4, "B")

# motif S4 (no. of apparent competition motifs/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s4, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.35,1), xlab="interval", ylab="S4: apparent competition")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 1, "C")

# motif S5 (no. of direct competition motifs/size)
par(mar=c(5,5,1,2))
plot(guilds$time, guilds$s5, type="l", lty=1, lwd=4, col=2, 
     axes=FALSE, ylim=c(0.15,0.35), xlab="interval", ylab="S5: direct competition")
axis(2)
axis(1, at = 1:4, labels =c("pre-extinction", "post-extinction", "early recovery", "late recovery"))
box()
text(4, 0.35, "D")

par(mfrow=c(1,1))