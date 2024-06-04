
import xml.etree.ElementTree as ET
import numpy as np


def get_dimension_value(dimension_str):
    """Convert SVG dimension (with 'px' or 'pt') to a float."""
    if 'px' in dimension_str:
        return float(dimension_str.replace('px', ''))
    elif 'pt' in dimension_str:
        # Assuming 1pt is approximately 1.33px (default for browsers).
        return float(dimension_str.replace('pt', '')) * 1.33
    else:
        return float(dimension_str)


def combine_svgs(svg_paths, shape, order, label, output_path):
    # Validate input
    if len(svg_paths) > shape[0] * shape[1]:
        print(
            f"There are only {shape[0] * shape[1]} spots for shape {shape}, but got {len(svg_paths)} svg files.")
        return

    # Map A,B,C.. to Number
    ABC = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
           'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    # Number = np.arange(1, shape[0]*shape[1]+1, 1)
    # ABC_to_Number = {}
    # for abc, num in zip(ABC, Number):
    #     dict = {abc: num}
    #     ABC_to_Number.update(dict)
    # # print(ABC_to_Number)

    # num_list = [ABC_to_Number[abc.split('.')[0]] for abc in order]

    num_list = order

    i = -1
    widths = []
    heights = []
    roots = []
    trees = []
    for _ in np.arange(1, shape[0]*shape[1]+1, 1):
        if _ in num_list:
            i += 1
            svg_path = svg_paths[i]
            tree = ET.parse(svg_path)
            root = tree.getroot()
            trees.append(tree)
            roots.append(root)
            widths.append(get_dimension_value(root.get('width')))
            heights.append(get_dimension_value(root.get('height')))
        else:
            trees.append('None')
            roots.append('None')
            widths.append(get_dimension_value(root.get('width')))
            heights.append(get_dimension_value(root.get('height')))

    total_width = [sum(widths[_*shape[1]:(_+1)*shape[1]])
                   for _ in range(shape[0])]
    total_width = np.max(total_width)

    total_height = [sum(heights[_*shape[0]:(_+1)*shape[0]])
                    for _ in range(shape[1])]
    total_height = np.max(total_height)

    # Create the new root SVG element
    new_root = ET.Element('svg', xmlns="http://www.w3.org/2000/svg", version="1.1",
                          width=str(total_width) + 'px', height=str(total_height) + 'px')

    # Position and append SVG elements to the new root
    _ = -1
    for i in range(shape[0]):
        for j in range(shape[1]):
            index = i * shape[1] + j
            root = roots[index]
            if root != 'None':
                _ += 1
                x_offset = sum(widths[k]-80 for k in range(j))
                y_offset = sum(heights[k * shape[1]]-45 for k in range(i))
                texts = label[_].split('//')
                if len(texts) == 2:
                    label_attributes = {
                        "x": str(x_offset+25),
                        "y": str(y_offset+40),
                        "font-family": "Arial",
                        "font-size": "40", }
                    label_ = ET.Element(
                        "{http://www.w3.org/2000/svg}text", **label_attributes)
                    label_.text = f"{texts[0]}"  # You can customize this text
                    new_root.append(label_)

                    label_attributes = {
                        "x": str(x_offset+widths[j]/2-40),
                        "y": str(y_offset+40),
                        "font-family": "Arial",
                        "font-size": "40", }
                    label_ = ET.Element(
                        "{http://www.w3.org/2000/svg}text", **label_attributes)
                    label_.text = f"{texts[1]}"  # You can customize this text
                    new_root.append(label_)

                if len(texts) == 1:
                    label_attributes = {
                        "x": str(x_offset+25),
                        "y": str(y_offset+40),
                        "font-family": "Arial",
                        "font-size": "40", }
                    label_ = ET.Element(
                        "{http://www.w3.org/2000/svg}text", **label_attributes)
                    label_.text = f"{texts[0]}"  # You can customize this text
                    new_root.append(label_)

                for element in root:
                    if element.tag.endswith(('g', 'path', 'circle', 'rect')):
                        if 'transform' in element.attrib:
                            transformations = element.attrib['transform'].split(
                                ' ')
                            for m, transform in enumerate(transformations):
                                if 'translate' in transform:
                                    values = transform.strip(
                                        'translate()').split(',')
                                    values[0] = str(
                                        float(values[0]) + x_offset+30)
                                    values[1] = str(
                                        float(values[1]) + y_offset+45)
                                    transformations[m] = 'translate(' + \
                                        ','.join(values) + ')'
                            element.attrib['transform'] = ' '.join(
                                transformations)
                        else:
                            element.attrib['transform'] = f"translate({x_offset+30},{y_offset+45})"
                    new_root.append(element)

    # Write the total SVG to the output path
    tree = ET.ElementTree(new_root)
    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    print('Finished.\n')


if __name__ == "__main__":
    import sys

    try:
        shape = tuple(map(int, sys.argv[1].split(',')))
        order = sys.argv[2]
        output = sys.argv[3]
        svg_paths = sys.argv[4:]
        print('Will array %d svgs with a shape of %s and save to %s' %
              (len(svg_paths), str(shape), output))

        if len(svg_paths) != shape[0] * shape[1]:
            raise ValueError
    except ValueError:
        print(
            f"There are only {shape[0] * shape[1]} spots for shape {shape}, but got {len(svg_paths)} svg files.")
        sys.exit()
    except:
        print('Usage: python combine_svgs.py <shape> <order> <output_path> <svg_path1> <svg_path2> ...')
        sys.exit()

    combine_svgs(svg_paths, shape, order, output)
