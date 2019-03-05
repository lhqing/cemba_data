def generate_tree_dict(data):
    """
    Helper function for echarts tree structures

    Parameters
    ----------
    data

    Returns
    -------

    """
    childrens = []
    if data.shape[1] == 1:
        cell_counts = data.iloc[:, 0].value_counts()
        for children, count in cell_counts.items():
            childrens.append({
                'name': children,
                'value': count
            })
    else:
        for children, sub_df in data.groupby(data.columns[0]):
            count = sub_df.shape[0]
            childrens.append({
                'name': children,
                'value': count,
                'children': generate_tree_dict(sub_df.iloc[:, 1:])
            })
    return childrens
