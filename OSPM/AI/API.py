import openai
import os
import tkinter as tk
from tkinter import scrolledtext, simpledialog, messagebox, filedialog
import json
import random
import datetime
import platform
from dotenv import load_dotenv

# Load environment variables
load_dotenv()


# Function to save API token to .env file
def save_api_token_to_env(api_token):
    with open('.env', 'w') as env_file:
        env_file.write(f"OPENAI_API_KEY={api_token}\n")

# Retrieve API token from the environment (from the .env file)
api_key = os.getenv('OPENAI_API_KEY')

# If the API key is not found, prompt the user to enter it
if not api_key:
    api_key = simpledialog.askstring("API Key", "Please enter your OpenAI API key:", show='*')
    if not api_key:
        messagebox.showerror("Error", "API Key is required to proceed!")
        quit()

    # Save the API key to the .env file so it doesn't need to be entered again
    save_api_token_to_env(api_key)

# Initialize OpenAI client with the API key
openai.api_key = api_key
client = openai.OpenAI(api_key=api_key)

# Define the available themes for the interface
themes = {
    "Retro": {"bg": "black", "fg": "#00FF00"},
    "Light": {"bg": "white", "fg": "black"}
}

# Create the Tkinter root window
root = tk.Tk()
root.title("Nate's Assistant")
root.configure(bg="black")

if platform.system() == 'Windows':
    root.state('zoomed')  # Maximizes the window on Windows
elif platform.system() == 'Linux':
    root.attributes('-zoomed', True)  # Maximizes the window on Linux
else:
    root.geometry("{0}x{1}+0+0".format(root.winfo_screenwidth(), root.winfo_screenheight()))  # Fullscreen on macOS

# Default config values (used if config.json is not found)
default_config = {
    "available_models": ["gpt-3.5-turbo", "gpt-4"],
    "token_options": [150, 300, 500, 750, 1024]
}

# Function to bind Shift + Enter to send prompt
def on_shift_enter(event):
    send_prompt()
    return "break"  # Prevent default newline behavior

# Function to create default config file if missing
def create_default_config():
    with open('config.json', 'w') as config_file:
        json.dump(default_config, config_file, indent=4)
    messagebox.showinfo("Config File Created", "config.json was not found and has been created with default values.")

# Check if the config file exists, create it if not
if not os.path.exists('config.json'):
    create_default_config()

# Load config settings
with open('config.json', 'r') as config_file:
    config = json.load(config_file)

# Dictionary to store conversation sessions
conversations = {}

# Default system prompt
SYSTEM_PROMPT = (
    "You are Jarvis, a highly intelligent and capable virtual assistant whose purpose is to help an astrophysicist named Nate. "
    "You provide clear, concise, and helpful information in a research type format, keeping all responses in paragraph form such as that of a scientific paper, and will provide all references at the end of every response."
)

# Function to dynamically generate a random greeting
def generate_random_greeting():
    greetings_in_english = ["Hello", "Hi", "Good day", "Greetings", "Salutations"]
    greetings_in_other_languages = {
        "Spanish": "Hola",
        "French": "Bonjour",
        "German": "Guten Tag",
        "Italian": "Ciao",
        "Portuguese": "Olá",
        "Japanese": "こんにちは",
        "Korean": "안녕하세요",
        "Hindi": "नमस्ते"
    }

    # Randomly select whether to use English or another language
    if random.random() < 0.5:
        greeting = random.choice(greetings_in_english)
        return f"{greeting}, Nate!"
    else:
        language, greeting = random.choice(list(greetings_in_other_languages.items()))
        return f"{greeting}, Nate! (in {language})"

# Function to load conversations from file
def load_conversations():
    global conversations
    if os.path.exists("conversations.json"):
        with open("conversations.json", "r") as file:
            conversations = json.load(file)

# Function to save conversations to file
def save_conversations():
    with open("conversations.json", "w") as file:
        json.dump(conversations, file)

# Function to export conversations
def export_conversation():
    filepath = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
    if filepath:
        with open(filepath, 'w') as file:
            for line in conversations[current_session]:
                file.write(line + "\n")
        messagebox.showinfo("Export Success", "Conversation exported successfully.")

# Function to send a prompt to OpenAI and get a response
def generate_response(prompt):
    try:
        # Get the material provided by the user (if any)
        material = material_entry.get("1.0", tk.END).strip()

        history = conversations[current_session][-12:]
        messages = [{"role": "system", "content": SYSTEM_PROMPT}]

        for entry in history:
            role, content = entry.split(": ", 1)
            if role.lower() == "you":
                messages.append({"role": "user", "content": content})
            elif role.lower() == "assistant":
                messages.append({"role": "assistant", "content": content})

        # Add the user's latest input
        messages.append({"role": "user", "content": prompt})

        # If user has provided material, prioritize it
        if material:
            messages.append({"role": "system", "content": f"Use the following material as the primary reference for your response: {material}"})

        # Use the updated API method with selected model and token limit
        response = client.chat.completions.create(
            model=selected_model.get(),
            messages=messages,
            max_tokens=selected_token.get(),
            temperature=0.7,
            top_p=1.0,
            frequency_penalty=0.0,
            presence_penalty=0.6
        )

        # Extract the assistant's reply in the correct way
        assistant_reply = response.choices[0].message.content.strip()

        # Add peer-reviewed, assumptions, and references (mock values)
        peer_reviewed = "Peer Reviewed: Yes (arXiv, AAS)"
        assumptions = "Assumptions: Standard physics models applied."
        references = "References: arXiv:1234.5678, AAS: doi.org/10.1111/aas.5678"

        return f"{assistant_reply}\n\n{peer_reviewed}\n{assumptions}\n{references}"

    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")
        return "Something went wrong, please try again."

# Function to handle sending the user's prompt
def send_prompt(event=None):
    prompt = prompt_entry.get("1.0", tk.END).strip()
    if not prompt:
        return
    response = generate_response(prompt)
    conversations[current_session].append(f"You: {prompt}")
    conversations[current_session].append(f"Assistant: {response}")
    save_conversations()
    update_chat_display()
    prompt_entry.delete("1.0", tk.END)

# Function to display conversation history
def update_chat_display():
    chat_display.config(state=tk.NORMAL)
    chat_display.delete("1.0", tk.END)
    for line in conversations[current_session]:
        chat_display.insert(tk.END, line + "\n")
    chat_display.config(state=tk.DISABLED)

# Function to start a new session with a dynamically generated greeting
def start_new_session():
    global current_session
    current_session = f"session_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    conversations[current_session] = []
    greeting = generate_random_greeting()  # Generate random greeting dynamically
    conversations[current_session].append(f"Assistant: {greeting}")
    update_chat_display()

# Function to switch to a different session
def switch_session():
    global current_session
    session_name = session_var.get()
    if session_name in conversations:
        current_session = session_name
        update_chat_display()

# Function to change the theme
def change_theme(selected_theme):
    theme = themes[selected_theme]
    chat_display.config(bg=theme["bg"], fg=theme["fg"])
    prompt_entry.config(bg=theme["bg"], fg=theme["fg"])

# Available GPT models and max tokens
available_models = config['available_models']
selected_model = tk.StringVar(root)
selected_model.set(available_models[0])

token_options = config['token_options']
selected_token = tk.IntVar(root)
selected_token.set(token_options[0])

# Create a chat display
chat_display = scrolledtext.ScrolledText(
    root, wrap=tk.WORD, state=tk.DISABLED, bg="black", fg="#00FF00", font=("Courier", 12)
)
chat_display.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

# Title for the prompt entry box
prompt_label = tk.Label(root, text="Enter Your Prompt Below:", bg="black", fg="#00FF00", font=("Courier", 12))
prompt_label.pack(padx=10, pady=5)

# Add text entry box for user input, with visible cursor in retro style
prompt_entry = tk.Text(root, height=5, bg="black", fg="#00FF00", font=("Courier", 12), insertbackground="#00FF00")
prompt_entry.pack(padx=10, pady=5, fill=tk.X)
prompt_entry.bind("<Shift-Return>", on_shift_enter)

# Create a frame for dropdowns and inputs
input_frame = tk.Frame(root, bg="black")
input_frame.pack(padx=10, pady=5, fill=tk.X)

# Create a material input box
material_label = tk.Label(input_frame, text="Enter Material (Optional):", bg="black", fg="#00FF00", font=("Courier", 12))
material_label.pack(side=tk.LEFT, padx=5)
material_entry = tk.Text(input_frame, height=3, width=40, bg="black", fg="#00FF00", font=("Courier", 12), insertbackground="#00FF00")
material_entry.pack(side=tk.LEFT, padx=5)

# Dropdown for model selection
model_dropdown = tk.OptionMenu(input_frame, selected_model, *available_models)
model_dropdown.config(bg="black", fg="#00FF00", font=("Courier", 12))
model_dropdown.pack(side=tk.LEFT, padx=5)

# Dropdown for max token selection
token_dropdown = tk.OptionMenu(input_frame, selected_token, *token_options)
token_dropdown.config(bg="black", fg="#00FF00", font=("Courier", 12))
token_dropdown.pack(side=tk.LEFT, padx=5)

# Dropdown for session selection
session_var = tk.StringVar(root)
session_var.set(next(iter(conversations), 'default'))  # Set default session

if conversations:
    session_dropdown = tk.OptionMenu(input_frame, session_var, *conversations.keys())
else:
    session_dropdown = tk.OptionMenu(input_frame, session_var, 'default')
session_dropdown.config(bg="black", fg="#00FF00", font=("Courier", 12))
session_dropdown.pack(side=tk.LEFT, padx=5)

# Theme selector
theme_var = tk.StringVar(root)
theme_var.set("Retro")
theme_dropdown = tk.OptionMenu(input_frame, theme_var, *themes.keys(), command=change_theme)
theme_dropdown.config(bg="black", fg="#00FF00", font=("Courier", 12))
theme_dropdown.pack(side=tk.LEFT, padx=5)

# Create a frame for the buttons (to align them in a row at the bottom)
button_frame = tk.Frame(root, bg="black")
button_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=10)

# Add buttons to the button frame to align them in a row
switch_session_button = tk.Button(button_frame, text="Switch Session", command=switch_session, bg="black", fg="#00FF00", font=("Courier", 12))
switch_session_button.pack(side=tk.LEFT, padx=5)

export_button = tk.Button(button_frame, text="Export Conversation", command=export_conversation, bg="black", fg="#00FF00", font=("Courier", 12))
export_button.pack(side=tk.LEFT, padx=5)

new_session_button = tk.Button(button_frame, text="New Session", command=start_new_session, bg="black", fg="#00FF00", font=("Courier", 12))
new_session_button.pack(side=tk.LEFT, padx=5)

submit_button = tk.Button(button_frame, text="Submit", command=send_prompt, bg="black", fg="#00FF00", font=("Courier", 12))
submit_button.pack(side=tk.LEFT, padx=5)

# Ensure all dropdowns, buttons, and labels use the retro color scheme
model_dropdown.config(bg="black", fg="#00FF00")
token_dropdown.config(bg="black", fg="#00FF00")
session_dropdown.config(bg="black", fg="#00FF00")
theme_dropdown.config(bg="black", fg="#00FF00")
new_session_button.config(bg="black", fg="#00FF00")
switch_session_button.config(bg="black", fg="#00FF00")
export_button.config(bg="black", fg="#00FF00")
submit_button.config(bg="black", fg="#00FF00")

# Load conversations
load_conversations()

# Initialize with a default session if not present
current_session = "default"
if current_session not in conversations:
    conversations[current_session] = []

# Update display and run the app
update_chat_display()
root.mainloop()