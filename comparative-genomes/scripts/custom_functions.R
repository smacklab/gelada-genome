#!/usr/bin/env Rscript

# Function to find the target node of a given taxon (CAFE)
find.target.node = function(x) {
	node.id.out = gsub(paste0('^.+?',x,'<([0-9]+)>.*'),'\\1',node.ids)
	node.id.out
}

# Function to find the triple (target node, parent node, and sister node) (CAFE)
find.triple = function(x) {
	suppressMessages(require(stringr))

	node.integers.reduced = node.integers
	node.position = unlist(str_locate_all(node.integers.reduced,paste0('<',x,'>')))

	if (str_sub(node.integers.reduced,start=node.position[1]-1,end=node.position[1]-1) == ')') {
		# If the preceding character is a close parenthesis, then the node itself is an ancestral node and should be collapsed
		sub.string = str_sub(node.integers.reduced,end=node.position[1]-1)
		start.parentheses = unique(unlist(str_locate_all(sub.string,'\\(')))
		end.parentheses = unique(unlist(str_locate_all(sub.string,'\\)')))
		all.parentheses = sort(c(start.parentheses,end.parentheses),decreasing=TRUE)
		parentheses.sum = numeric(length(all.parentheses))
		parentheses.sum[match(start.parentheses,all.parentheses)] = -1
		parentheses.sum[match(end.parentheses,all.parentheses)] = 1
		node.integers.reduced = paste0(str_sub(sub.string,end=all.parentheses[which(!cumsum(parentheses.sum))[1]]-1),str_sub(node.integers.reduced,start=node.position[1]))
	}
	node.position = unlist(str_locate_all(node.integers.reduced,paste0('<',x,'>')))
	if (str_sub(node.integers.reduced,start=node.position[2]+1,end=node.position[2]+2) == ',(') {
		# If the following characters are a comma followed by an open parenthesis, then the sister is an ancestral node and should be collapsed
		sub.string = str_sub(node.integers.reduced,start=node.position[2]+2)
		start.parentheses = unique(unlist(str_locate_all(sub.string,'\\(')))
		end.parentheses = unique(unlist(str_locate_all(sub.string,'\\)')))
		all.parentheses = sort(c(start.parentheses,end.parentheses),decreasing=FALSE)
		parentheses.sum = numeric(length(all.parentheses))
		parentheses.sum[match(start.parentheses,all.parentheses)] = 1
		parentheses.sum[match(end.parentheses,all.parentheses)] = -1
		node.integers.reduced = paste0(str_sub(node.integers.reduced,end=node.position[2]+1),str_sub(sub.string,start=all.parentheses[which(!cumsum(parentheses.sum))[1]]+1))
	} else if (length(unlist(str_locate_all(node.integers.reduced,paste0('\\)<[0-9]+>,<',x,'>'))))) {
		# If the preceding characters are a closed parenthesis followed a node ID and a comma, then the sister is an ancestral node and should be collapsed
		node.positions = unlist(str_locate_all(node.integers.reduced,paste0('\\)<[0-9]+>,<',x,'>')))
		sub.string = str_sub(node.integers.reduced,end=node.positions[1])
		start.parentheses = unique(unlist(str_locate_all(sub.string,'\\(')))
		end.parentheses = unique(unlist(str_locate_all(sub.string,'\\)')))
		all.parentheses = sort(c(start.parentheses,end.parentheses),decreasing=TRUE)
		parentheses.sum = numeric(length(all.parentheses))
		parentheses.sum[match(start.parentheses,all.parentheses)] = -1
		parentheses.sum[match(end.parentheses,all.parentheses)] = 1
		node.integers.reduced = paste0(str_sub(sub.string,end=all.parentheses[which(!cumsum(parentheses.sum))[1]]-1),str_sub(node.integers.reduced,start=node.positions[1]+1))
	}

	triple = str_extract(node.integers.reduced,paste0('\\([0-9<>,]*<',x,'>[0-9<>,]*\\)<[0-9]+>'))
	if (is.na(triple)) {
		str_extract(node.integers.reduced,paste0('<',x,'>'))
	} else {
		# Ensure that the target is shown first
		gsub(paste0('\\((.*?),(<',x,'>)\\)'),'(\\2,\\1)',triple)
	}
}

# Function to find the parent node of a given taxon (CAFE)
find.parent.node = function(x) {
	suppressMessages(require(stringr))

	triple = find.triple(x)
	
	if (str_count(triple,'<[0-9]+>') < 3) {
		NA
	} else {
		gsub('\\(<([0-9]+)>,<([0-9]+)>\\)<([0-9]+)>','\\3',triple)
	}
}

# Function to find the sister node of a given taxon (CAFE)
find.sister.node = function(x) {
	suppressMessages(require(stringr))

	triple = find.triple(x)
	
	if (str_count(triple,'<[0-9]+>') < 3) {
		NA
	} else {
		gsub('\\(<([0-9]+)>,<([0-9]+)>\\)<([0-9]+)>','\\2',triple)
	}
}

# Function to find the index of the p-value given a node ID (CAFE)
find.node.index = function(x) {
	match(x,unlist(strsplit(gsub(';$','',gsub('^;','',gsub('[(), ]+',';',unlist(strsplit(scan(results.file,quiet=TRUE,nlines=1,skip=3,what='',sep='\n'),':'))[3]))),';')))
}

# Function to find the protein count of a given node

get.node.count = function(tree,x) {
	as.integer(gsub('^([0-9]+).*','\\1',gsub(gsub('([()])','\\\\\\1',unlist(strsplit(gsub('<[0-9]+>','_[0-9]+',gsub(paste0('<',x,'>'),'_#',node.ids)),'#'))[1]),'',gsub(':[0-9.]+','',tree))))
}

# Function to find the p-value of a given node (CAFE)

get.node.p = function(p,i) {
	out = unlist(strsplit(p,'[(),]+'))
	out = as.numeric(out[as.logical(nchar(out))])
	out[i]
}

# Function to parse a newick and remove a taxon (adjusting branch lengths)

prune.branch = function(x,nwk) {
	suppressMessages(require(stringr))
	branch.position = unlist(str_locate_all(nwk,paste0(x,':[0-9.]+')))
	branch = str_sub(nwk,start=matrix(branch.position,nrow=1))
	branch.length = as.numeric(gsub('.+?:([0-9.]+)$','\\1',branch))

	# If preceding comma
	if (str_sub(nwk,start=branch.position[1]-1,end=branch.position[1]-1) == ',') {
		sub.string = str_sub(nwk,end=branch.position[1]-2)
		start.parentheses = unique(unlist(str_locate_all(sub.string,'\\(')))
		end.parentheses = unique(unlist(str_locate_all(sub.string,'\\)')))
		all.parentheses = sort(c(start.parentheses,end.parentheses),decreasing=TRUE)
		parentheses.sum = numeric(length(all.parentheses))
		parentheses.sum[match(start.parentheses,all.parentheses)] = 1
		parentheses.sum[match(end.parentheses,all.parentheses)] = -1
		sister = str_sub(sub.string,start=all.parentheses[which(cumsum(parentheses.sum) == 1)[1]]+1)
		node.positions = unlist(str_locate_all(nwk,gsub('([()])','\\\\\\1',paste0('(',sister,',',branch,')',':[0-9.]+'))))
	} else if (str_sub(nwk,start=branch.position[2]+1,end=branch.position[2]+1) == ',') {
		sub.string = str_sub(nwk,start=branch.position[2]+2)
		start.parentheses = unique(unlist(str_locate_all(sub.string,'\\(')))
		end.parentheses = unique(unlist(str_locate_all(sub.string,'\\)')))
		all.parentheses = sort(c(start.parentheses,end.parentheses),decreasing=FALSE)
		parentheses.sum = numeric(length(all.parentheses))
		parentheses.sum[match(start.parentheses,all.parentheses)] = -1
		parentheses.sum[match(end.parentheses,all.parentheses)] = 1
		sister = str_sub(sub.string,end=all.parentheses[which(cumsum(parentheses.sum) == 1)[1]]-1)
		node.positions = unlist(str_locate_all(nwk,gsub('([()])','\\\\\\1',paste0('(',branch,',',sister,')',':[0-9.]+'))))
	}

	if (length(node.positions)) {
		node = str_sub(nwk,matrix(node.positions,nrow=1))
		sister.branch.length = as.numeric(gsub('.+?:([0-9.]+)$','\\1',sister))
		root.length = as.numeric(gsub('.+?:([0-9.]+)$','\\1',node))
		new.node = gsub('(.+?:)[0-9.]+$',paste0('\\1',sister.branch.length+root.length),sister)
		nwk.out = paste0(str_sub(nwk,end=node.positions[1]-1),new.node,str_sub(nwk,start=node.positions[2]+1))
	} else {
		nwk.out = gsub(':[0-9.]+$',';',sister)
	}
	nwk.out
}