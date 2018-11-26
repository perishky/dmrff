region.members <- function(intervals, positions) {
    stopifnot(is.data.frame(intervals))
    stopifnot(all(c("chr","start","end") %in% colnames(intervals)))
    stopifnot(is.data.frame(positions))
    stopifnot(all(c("chr","pos") %in% colnames(positions)))
    
    events <- rbind(data.frame(chr=as.character(intervals$chr),
                               pos=intervals$start,
                               type="start",
                               id=1:nrow(intervals)),
                    data.frame(chr=as.character(intervals$chr),
                               pos=intervals$end,
                               type="end",
                               id=1:nrow(intervals)),
                    data.frame(chr=as.character(positions$chr),
                               pos=positions$pos,
                               type="pos",
                               id=1:nrow(positions)))
    events <- events[order(events$chr, events$pos, sign(events$type!="start"), sign(events$type=="end"), decreasing=F),]
    start.idx <- which(events$type == "start")
    end.idx <- which(events$type == "end")
    intervals$event.start.idx <- start.idx[match(1:nrow(intervals), events$id[start.idx])] + 1
    intervals$event.end.idx <- end.idx[match(1:nrow(intervals), events$id[end.idx])] - 1
    events$id[which(events$type != "pos")] <- NA    
    lapply(1:nrow(intervals), function(i) {
        if (intervals$event.start.idx[i] >= intervals$event.end.idx[i]) integer(0)
        na.omit(events$id[intervals$event.start.idx[i]:intervals$event.end.idx[i]])
    })         
}
