"""
# Description: Contains utility functions for beautifying string outputs.

# Functions covered:
    - print_msg_box: Prints a message box around a given message with an optional title.

# Dependencies:
    - None
"""

def print_msg_box(msg, indent=1, width=None, title=None):
    """
    Prints a message box around a given message with an optional title.

    Args:
        - msg (str): The message to display inside the box. It can contain newlines.
        - indent (int): The number of spaces to indent the box contents. Default is 1.
        - width (int, optional): The width of the message box. If not provided, the width is automatically set to the 
                            length of the longest line in the message.
        - title (str, optional): A title to be displayed at the top of the message box. If provided, it will appear above the message.

    Returns:
        - None: This function only prints the message box to the console, it does not return any value.
    """
    
    # Split the message into individual lines based on newlines
    lines = msg.split('\n')
    
    # Create a string representing the indentation (spaces before content)
    space = " " * indent
    
    # If no width is provided, set it to the length of the longest line in the message
    if not width:
        width = max(map(len, lines))
    
    # Initialize the top border of the box, including the indentation
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    
    # If a title is provided, add the title section to the box
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore under the title
    
    # Add the message lines to the box, ensuring they are left-aligned within the box width
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    
    # Add the bottom border of the box
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    
    # Print the final message box to the console
    print(box)
