"""Display names only; original annotation keys remain unchanged for scoring."""
DISPLAY_NAMES = {
    'P-TGC (Sec. P-TGC 1)': 'Secondary P-TGC1',
    'P-TGC (Sec. P-TGC 2)': 'Secondary P-TGC2',
}

def display_name(name):
    return DISPLAY_NAMES.get(name,name)

def rename_reference(root):
    """Retain the two-line layout at the existing nodes, without parentheses."""
    changes=0
    for node in root.iter():
        if not node.tag.endswith('}text') or node.text is None:continue
        s=''.join(node.text.split())
        if s=='P-TGC':node.text='Secondary';changes+=1
        elif s in ['(Sec.P-TGC1)','(Sec.P-TGC2)']:
            node.text='P-TGC'+s[-2]
            node.set('font-size','5.8');changes+=1
    return changes
