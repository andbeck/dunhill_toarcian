library(cheddar)

data(TL84) # Tuesday Lake

TL_nodes <- TL84$nodes$node
primary <- BasalNodes(TL84)
out <- RemoveNodes(TL84, BasalNodes(TL84), method='cascade')
end_nodes <- out$nodes$node

# secondary are not primary or end.
!(TL_nodes %in% primary) # 31 secondary AND survivors

intermediate <- TL_nodes[!(TL_nodes %in% primary)]

sum(!(intermediate %in% end_nodes)) # 25 secondaries

secondary<- intermediate[!(intermediate %in% end_nodes)]
