from types import ListType


entity_levels=["A", "R", "C", "M", "S"]

def uniqueify(l):
	"Return unique items in list l."
	d={}
	for i in l:
		if not d.has_key(i):
			d[i]=None
	return d.keys()

def get_unique_parents(entity_list):
	l=[]
	for entity in entity_list:
		parent=entity.get_parent()
		l.append(parent)
	return uniqueify(l)

def unfold_entities(entity_list, target_level):
	if not target_level in entity_levels:
		raise Exception, "%s: Not an entity level." % target_level
	if type(entity_list)!=ListType:
		entity_list=[entity_list]
	level=entity_list[0].get_level()
	for entity in entity_list:
		if not (entity.get_level()==level):
			raise Exception, "Entity list is not homogeneous."
	target_index=entity_levels.index(target_level)
	level_index=entity_levels.index(level)
	if level_index==target_index:
		return [entity]
	if target_index==level_index-1:
		return entity.get_list()
	if target_index==level_index+1:
		return [entity.get_parent()]
	if level_index>target_index:
		for i in range(target_index, level_index):
			new_entity_list=[]
			for entity in entity_list:
				new_entity_list=new_entity_list+entity.get_list()
			entity_list=new_entity_list
	else:
		for i in range(level_index, target_index):
			new_entity_list=[]	
			for entity in entity_list:
				parent=entity.get_parent()
				new_entity_list.append(parent)
			entity_list=uniqueify(new_entity_list)
	return entity_list

