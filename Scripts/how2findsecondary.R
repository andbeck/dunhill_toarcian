library(cheddar)

data(TL84) # Tuesday Lake

start <- TL84$nodes$node
primary <- BasalNodes(TL84)
out <- RemoveNodes(TL84, BasalNodes(TL84), method='cascade')
end <- out$nodes$node

# secondary are not primary or end.
!(start %in% primary) # 31 secondary AND survivors

intermediate <- start[!(start %in% primary)]

sum(!(intermediate %in% end)) # 25 secondaries

secondary<- intermediate[!(intermediate %in% end)]


start[!(start %in% primary)][!(intermediate %in% end)]
