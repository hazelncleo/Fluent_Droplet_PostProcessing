from inquirer.themes import Default
from blessed import Terminal

term = Terminal()

def green_text(text):
    return term.green3(str(text))

def red_text(text):
    return term.firebrick1(str(text))

def blue_text(text):
    return term.deepskyblue(str(text))

def yellow_text(text):
    return term.gold(str(text))

class HazelsAwesomeTheme(Default):
    def __init__(self):
        super().__init__()
        self.Question.mark_color = term.deepskyblue
        self.Question.brackets_color = term.normal
        self.Question.default_color = term.normal
        self.Editor.opening_prompt_color = term.gray33
        self.Checkbox.selection_color = term.deepskyblue
        self.Checkbox.selection_icon = ">"
        self.Checkbox.selected_icon = "[X]"
        self.Checkbox.selected_color = term.gold + term.bold
        self.Checkbox.unselected_color = term.normal
        self.Checkbox.unselected_icon = "[ ]"
        self.Checkbox.locked_option_color = term.gray50
        self.List.selection_color = term.deepskyblue
        self.List.selection_cursor = ">"
        self.List.unselected_color = term.normal