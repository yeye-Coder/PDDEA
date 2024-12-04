import tkinter as tk
from tkinter import filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Rectangle
import os

class PDDEA:
    def __init__(self, master):
        self.master = master
        self.master.title("Partial Discharge Dataset Editing Application (PDDEA)")

        # Initialize variables
        self.prpd = None
        self.waveform = None
        self.peak_max = None
        self.current_sn = 0  # Serial Number
        self.selected_dot = None
        self.group_selected_indices = []
        self.group_selected_dots = []

        # Variables for rectangle selection
        self.press_event = None
        self.dragging = False
        self.drag_threshold = 5  # pixels, adjust as necessary
        self.rect_select = None

        # Variables to store original filenames
        self.prpd_filename = ''
        self.waveform_filename = ''

        # Variable for the "Negate Waveform" checkbox
        self.negate_waveform_var = tk.BooleanVar(value=False)

        # Configure the master grid
        self.master.columnconfigure(0, weight=1)
        self.master.columnconfigure(1, weight=1)

        # Top row components
        self.create_top_row()

        # Second row components (Graph panes)
        self.create_graph_panes()

        # Create context menu
        self.create_context_menu()

        # Add footer label
        self.create_footer_label()

    def create_top_row(self):
        top_frame = tk.Frame(self.master)
        top_frame.grid(row=0, column=0, columnspan=2, sticky='ew')
        top_frame.columnconfigure(0, weight=1)
        top_frame.columnconfigure(1, weight=1)

        # Left subframe for labels and checkbox
        label_frame = tk.Frame(top_frame)
        label_frame.grid(row=0, column=0, sticky='w', padx=10)

        # Right subframe for buttons
        button_frame = tk.Frame(top_frame)
        button_frame.grid(row=0, column=1, sticky='e', padx=10)

        # Labels
        tk.Label(label_frame, text="X:").pack(side=tk.LEFT, padx=5)
        self.x_label = tk.Label(label_frame, text="0.0000", width=8, anchor='w')
        self.x_label.pack(side=tk.LEFT, padx=5)

        tk.Label(label_frame, text="Y:").pack(side=tk.LEFT, padx=5)
        self.y_label = tk.Label(label_frame, text="0.0000", width=8, anchor='w')
        self.y_label.pack(side=tk.LEFT, padx=5)

        tk.Label(label_frame, text="SN:").pack(side=tk.LEFT, padx=5)
        self.sn_label = tk.Label(label_frame, text="0", width=5, anchor='w')
        self.sn_label.pack(side=tk.LEFT, padx=5)

        # "Negate Waveform" checkbox
        self.negate_waveform_checkbox = tk.Checkbutton(
            label_frame,
            text="Negate Waveform",
            variable=self.negate_waveform_var
        )
        self.negate_waveform_checkbox.pack(side=tk.LEFT, padx=10)

        # Buttons
        tk.Button(button_frame, text="Open", command=self.open_files).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Save", command=self.save_files).pack(side=tk.LEFT, padx=5)
        tk.Button(button_frame, text="Exit", command=self.master.quit).pack(side=tk.LEFT, padx=5)

    def create_graph_panes(self):
        self.canvas_frame = tk.Frame(self.master)
        self.canvas_frame.grid(row=1, column=0, columnspan=2, sticky='nsew')
        self.master.rowconfigure(1, weight=1)
        self.canvas_frame.columnconfigure(0, weight=1)
        self.canvas_frame.columnconfigure(1, weight=1)

        # Left Pane (PRPD)
        self.prpd_fig, self.prpd_ax = plt.subplots(figsize=(6, 4))
        self.prpd_canvas = FigureCanvasTkAgg(self.prpd_fig, master=self.canvas_frame)
        self.prpd_canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        self.prpd_canvas.mpl_connect("motion_notify_event", self.on_prpd_motion)
        self.prpd_canvas.mpl_connect("button_press_event", self.on_prpd_press)
        self.prpd_canvas.mpl_connect("button_release_event", self.on_prpd_release)

        # Right Pane (Waveform)
        self.waveform_fig, self.waveform_ax = plt.subplots(figsize=(6, 4))
        self.waveform_canvas = FigureCanvasTkAgg(self.waveform_fig, master=self.canvas_frame)
        self.waveform_canvas.get_tk_widget().grid(row=0, column=1, sticky='nsew')

    def create_footer_label(self):
        # Footer label with corrected text
        footer_label = tk.Label(
            self.master,
            text="© LHQ & ChatGPT o1 2024",
            font=("Arial", 10),
            anchor='center'
        )
        footer_label.grid(row=2, column=0, columnspan=2, pady=5)

    def create_context_menu(self):
        # Create a context menu
        self.context_menu = tk.Menu(self.master, tearoff=0)
        self.context_menu.add_command(label="Delete", command=self.delete_points)
        self.context_menu.add_command(label="Negate", command=self.negate_peaks)

    def open_files(self):
        while True:
            # Select PRPD and Waveform files
            files = filedialog.askopenfilenames(
                title="Select PRPD and Waveform Files",
                filetypes=[("CSV files", "*.csv")]
            )
            if not files:
                return  # User canceled the dialog

            if len(files) != 2:
                messagebox.showerror("Error", "Please select exactly two files: PRPD and Waveform CSV files.")
                continue  # Prompt again

            # Initialize variables
            prpd_file = None
            waveform_file = None
            filenames_valid = True

            for file in files:
                filename = os.path.basename(file)
                filename_lower = filename.lower()
                if 'prpd' in filename_lower:
                    prpd_file = file
                    self.prpd_filename = filename
                elif 'waveform' in filename_lower:
                    waveform_file = file
                    self.waveform_filename = filename
                else:
                    messagebox.showerror("Error", "Each filename must contain 'prpd' or 'waveform'.")
                    filenames_valid = False
                    break  # Exit the for loop

            if not filenames_valid:
                continue  # Prompt again

            if not prpd_file or not waveform_file:
                messagebox.showerror("Error", "Please select one 'prpd' file and one 'waveform' file.")
                continue  # Prompt again

            # Define a function to extract the base name
            def get_base_name(filename):
                # Remove 'prpd' and 'waveform' (case-insensitive), and remove underscores and hyphens
                filename_lower = filename.lower()
                base = filename_lower.replace('prpd', '').replace('waveform', '')
                # Remove underscores, hyphens, and spaces
                base = base.replace('_', '').replace('-', '').replace(' ', '')
                # Remove the file extension
                base = os.path.splitext(base)[0]
                return base

            # Get base names
            base_name_prpd = get_base_name(self.prpd_filename)
            base_name_waveform = get_base_name(self.waveform_filename)

            if base_name_prpd != base_name_waveform:
                messagebox.showerror(
                    "Error",
                    "The filenames must match except for 'prpd' and 'waveform'.\n"
                    f"Base name mismatch: '{base_name_prpd}' vs '{base_name_waveform}'."
                )
                continue  # Prompt again

            try:
                # Helper function to load and clean data
                def load_and_clean(file_path, expected_rows=None):
                    with open(file_path, 'r') as file:
                        # Read and clean the data
                        data = file.read()
                        # Replace tabs and spaces with commas, remove newlines
                        data = data.replace('\n', ',').replace('\t', ',').replace(' ', '').strip()
                        # Split into numpy array
                        data = np.array([float(value) for value in data.split(',') if value.strip()], dtype=float)
                        # Reshape if required
                        if expected_rows:
                            data = data.reshape(expected_rows, -1)
                        return data

                # Load PRPD and waveform data
                self.prpd = load_and_clean(prpd_file, expected_rows=2)
                self.waveform = load_and_clean(waveform_file)

                # Validate PRPD and waveform shape relationships
                if self.waveform.shape[0] == self.prpd.shape[1] * 256:
                    self.waveform = self.waveform.reshape(self.prpd.shape[1], 256)
                else:
                    messagebox.showerror(
                        "Error",
                        f"Mismatch between PRPD and Waveform data.\n"
                        f"PRPD columns: {self.prpd.shape[1]}, Waveform rows: {self.waveform.shape[0]} (expected {self.prpd.shape[1] * 256})."
                    )
                    continue  # Prompt again

                # Debugging: Print shapes
                print(f"PRPD shape: {self.prpd.shape}")
                print(f"Waveform shape after reshaping: {self.waveform.shape}")

                # Validate shapes
                if self.prpd.shape[0] != 2:
                    messagebox.showerror("Error", "PRPD file must have 2 rows.")
                    continue  # Prompt again

                # Calculate peak_max
                self.peak_max = np.max(np.abs(self.waveform))
                if self.peak_max == 0:
                    self.peak_max = 1  # Avoid division by zero

                # Initialize SN
                self.current_sn = 0
                self.sn_label.config(text=str(self.current_sn))

                # Clear any previous selections
                self.group_selected_indices = []
                self.group_selected_dots = []

                # Plot PRPD and initial waveform
                self.plot_prpd()
                self.plot_waveform()
                break  # Successfully loaded, exit the loop

            except Exception as e:
                messagebox.showerror("Error", f"Failed to load files:\n{e}")
                continue  # Prompt again


    def plot_prpd(self):
        self.prpd_ax.clear()
        if self.prpd is not None and self.prpd.shape[1] > 0:
            self.peaks = self.prpd[0, :]
            self.phases = self.prpd[1, :]
            self.prpd_scatter = self.prpd_ax.scatter(self.phases, self.peaks, c='blue')
            self.prpd_ax.set_xlim(0, 360)
            self.prpd_ax.set_ylim(-self.peak_max, self.peak_max)
            self.prpd_ax.set_xlabel("Phase Angle (degrees)")
            self.prpd_ax.set_ylabel("Peak Value")
            self.prpd_ax.set_title("PRPD Pattern")

            # Highlight selected dot
            if self.selected_dot:
                self.selected_dot.remove()
                self.selected_dot = None

            if self.current_sn is not None and self.current_sn < len(self.phases):
                self.selected_dot, = self.prpd_ax.plot(
                    self.phases[self.current_sn],
                    self.peaks[self.current_sn],
                    'ro',
                    markersize=12,
                    fillstyle='none'
                )

            # Highlight group selected dots
            # Remove previous group selection circles
            for marker in getattr(self, 'group_selected_dots', []):
                marker.remove()
            self.group_selected_dots = []
            for idx in self.group_selected_indices:
                if idx < len(self.phases):
                    marker, = self.prpd_ax.plot(
                        self.phases[idx],
                        self.peaks[idx],
                        'go',
                        markersize=12,
                        fillstyle='none'
                    )
                    self.group_selected_dots.append(marker)

            # Remove rectangle selector if exists
            if self.rect_select:
                self.rect_select.remove()
                self.rect_select = None
        else:
            self.prpd_ax.set_title("No PRPD data")
            self.prpd_ax.set_xlim(0, 360)
            self.prpd_ax.set_ylim(-1, 1)
            self.prpd_ax.set_xlabel("Phase Angle (degrees)")
            self.prpd_ax.set_ylabel("Peak Value")

        self.prpd_canvas.draw()

    def plot_waveform(self):
        self.waveform_ax.clear()
        if self.waveform is not None and self.waveform.shape[0] > 0 and self.current_sn is not None:
            waveform_data = self.waveform[self.current_sn, :]
            n = waveform_data.size
            x = np.arange(n)
            self.waveform_ax.plot(x, waveform_data)
            self.waveform_ax.set_xlim(0, n - 1)
            self.waveform_ax.set_ylim(-self.peak_max, self.peak_max)
            self.waveform_ax.set_xlabel("Sample Index")
            self.waveform_ax.set_ylabel("Amplitude")
            self.waveform_ax.set_title(f"Waveform SN: {self.current_sn}")
        else:
            self.waveform_ax.set_title("No waveform data")
            self.waveform_ax.set_xlim(0, 1)
            self.waveform_ax.set_ylim(-1, 1)
            self.waveform_ax.set_xlabel("Sample Index")
            self.waveform_ax.set_ylabel("Amplitude")

        self.waveform_canvas.draw()

    def on_prpd_motion(self, event):
        if event.inaxes == self.prpd_ax and self.press_event:
            dx = event.x - self.press_event.x
            dy = event.y - self.press_event.y
            distance = (dx**2 + dy**2)**0.5
            if distance > self.drag_threshold:
                self.dragging = True
                if not self.rect_select:
                    # Initialize rectangle
                    self.rect_select = Rectangle(
                        (self.rect_start[0], self.rect_start[1]),
                        0, 0,
                        linewidth=1,
                        edgecolor='black',
                        linestyle='dotted',
                        facecolor='none'
                    )
                    self.prpd_ax.add_patch(self.rect_select)
                x0, y0 = self.rect_start
                x1, y1 = event.xdata, event.ydata
                if x1 is None or y1 is None:
                    return
                # Update rectangle dimensions
                self.rect_select.set_width(x1 - x0)
                self.rect_select.set_height(y1 - y0)
                self.rect_select.set_xy((x0, y0))
                self.prpd_canvas.draw()
        elif event.inaxes == self.prpd_ax:
            x, y = event.xdata, event.ydata
            if x is not None and y is not None:
                self.x_label.config(text=f"{x:.4f}")
                self.y_label.config(text=f"{y:.4f}")

    def on_prpd_press(self, event):
        if event.inaxes == self.prpd_ax:
            if event.button == 1:  # Left mouse button
                self.press_event = event
                self.dragging = False
                self.rect_start = (event.xdata, event.ydata)
            elif event.button == 3:  # Right mouse button
                self.show_context_menu(event)

    def on_prpd_release(self, event):
        if event.inaxes == self.prpd_ax and self.press_event:
            if self.dragging:
                # Handle rectangle selection
                x0, y0 = self.rect_start
                x1, y1 = event.xdata, event.ydata
                if x1 is None or y1 is None:
                    # Remove rectangle if mouse leaves the axes
                    if self.rect_select:
                        self.rect_select.remove()
                        self.rect_select = None
                    self.prpd_canvas.draw()
                    return
                # Normalize rectangle coordinates
                xmin = min(x0, x1)
                xmax = max(x0, x1)
                ymin = min(y0, y1)
                ymax = max(y0, y1)

                # Find points within rectangle
                indices = np.where(
                    (self.phases >= xmin) & (self.phases <= xmax) &
                    (self.peaks >= ymin) & (self.peaks <= ymax)
                )[0]

                # Update group selection
                self.group_selected_indices = indices.tolist()
                if indices.size > 0:
                    self.current_sn = indices[0]
                    self.sn_label.config(text=str(self.current_sn))
                else:
                    # If no dots are selected, clear group selection and keep current_sn unchanged
                    self.group_selected_indices = []

                # Update plots
                self.plot_prpd()
                self.plot_waveform()
            else:
                # Handle single click selection
                self.select_dot(event)

            # Reset flags
            self.press_event = None
            self.dragging = False

            # Remove rectangle selector if exists
            if self.rect_select:
                self.rect_select.remove()
                self.rect_select = None
            self.prpd_canvas.draw()

    def select_dot(self, event):
        x_click = event.xdata
        y_click = event.ydata
        if x_click is None or y_click is None:
            return

        # Define thresholds
        x_threshold = 5  # degrees
        y_threshold = self.peak_max * 0.05  # 5% of peak_max

        # Find dots within thresholds
        x_diff = np.abs(self.phases - x_click)
        y_diff = np.abs(self.peaks - y_click)

        # Handle circular phase angle (0 and 360 are the same)
        x_diff = np.minimum(x_diff, 360 - x_diff)

        close_points = np.where((x_diff <= x_threshold) & (y_diff <= y_threshold))[0]

        if close_points.size == 0:
            # No dot is close enough; do nothing
            return

        # Select the closest dot among the close points
        distances = np.hypot(x_diff[close_points], y_diff[close_points])
        nearest_index = close_points[np.argmin(distances)]

        # Update SN
        self.current_sn = nearest_index
        self.sn_label.config(text=str(self.current_sn))

        # Remove group selection markers
        for marker in getattr(self, 'group_selected_dots', []):
            marker.remove()
        self.group_selected_dots = []

        # Clear any group selection
        self.group_selected_indices = []

        # Update selected dot
        if self.selected_dot:
            self.selected_dot.remove()
        self.selected_dot, = self.prpd_ax.plot(
            self.phases[nearest_index],
            self.peaks[nearest_index],
            'ro',
            markersize=12,
            fillstyle='none'
        )

        self.prpd_canvas.draw()

        # Update waveform plot
        self.plot_waveform()

    def show_context_menu(self, event):
        x_click = event.xdata
        y_click = event.ydata
        if x_click is None or y_click is None:
            return

        # Define thresholds
        x_threshold = 5  # degrees
        y_threshold = self.peak_max * 0.05  # 5% of peak_max

        if self.group_selected_indices:
            # Check if right-click is on any of the group-selected dots
            x_diff = np.abs(self.phases[self.group_selected_indices] - x_click)
            y_diff = np.abs(self.peaks[self.group_selected_indices] - y_click)

            # Handle circular phase angle
            x_diff = np.minimum(x_diff, 360 - x_diff)

            close_points = np.where((x_diff <= x_threshold) & (y_diff <= y_threshold))[0]
            if close_points.size == 0:
                return  # Right-click is not on any selected dot
            else:
                # Show the context menu
                try:
                    self.context_menu.tk_popup(event.guiEvent.x_root, event.guiEvent.y_root)
                finally:
                    self.context_menu.grab_release()
        else:
            # Check if right-click is on the focused dot
            x_diff = np.abs(self.phases[self.current_sn] - x_click)
            y_diff = np.abs(self.peaks[self.current_sn] - y_click)

            # Handle circular phase angle
            x_diff = min(x_diff, 360 - x_diff)

            if (x_diff <= x_threshold) and (y_diff <= y_threshold):
                # Show the context menu
                try:
                    self.context_menu.tk_popup(event.guiEvent.x_root, event.guiEvent.y_root)
                finally:
                    self.context_menu.grab_release()
            else:
                return  # Right-click is not on the focused dot

    def delete_points(self):
        if self.group_selected_indices:
            indices = sorted(self.group_selected_indices, reverse=True)
            # Check if current_sn is among the indices being deleted
            deleting_current_sn = self.current_sn in indices

            for index in indices:
                self.prpd = np.delete(self.prpd, index, axis=1)
                self.waveform = np.delete(self.waveform, index, axis=0)
                if index < self.current_sn:
                    self.current_sn -= 1

            # Clear group selection
            self.group_selected_indices = []
            self.group_selected_dots = []

            # Update current_sn
            if deleting_current_sn:
                if self.prpd.shape[1] > 0:
                    self.current_sn = 0
                else:
                    self.current_sn = None  # No data left
        else:
            index = self.current_sn
            self.prpd = np.delete(self.prpd, index, axis=1)
            self.waveform = np.delete(self.waveform, index, axis=0)

            # Adjust current_sn if necessary
            if self.prpd.shape[1] > 0:
                if self.current_sn >= self.prpd.shape[1]:
                    self.current_sn = self.prpd.shape[1] - 1
            else:
                self.current_sn = None  # No data left

        # Recalculate peak_max
        if self.waveform.size > 0:
            self.peak_max = np.max(np.abs(self.waveform))
            if self.peak_max == 0:
                self.peak_max = 1  # Avoid division by zero
        else:
            self.peak_max = 1  # Avoid division by zero

        # Update sn_label
        if self.current_sn is not None:
            self.sn_label.config(text=str(self.current_sn))
        else:
            self.sn_label.config(text="N/A")

        # Update plots
        self.plot_prpd()
        self.plot_waveform()

    def negate_peaks(self):
        negate_waveform = self.negate_waveform_var.get()
        if self.group_selected_indices:
            indices = self.group_selected_indices
            # Negate the peak values
            self.prpd[0, indices] *= -1
            if negate_waveform:
                # Negate the waveform data
                self.waveform[indices, :] *= -1
            # Clear group selection
            self.group_selected_indices = []
            self.group_selected_dots = []
        else:
            index = self.current_sn
            # Negate the peak value
            self.prpd[0, index] *= -1
            if negate_waveform:
                # Negate the waveform data
                self.waveform[index, :] *= -1

        # Update peaks array for plotting
        self.peaks = self.prpd[0, :]

        # Update plots
        self.plot_prpd()
        self.plot_waveform()

    def save_files(self):
        if self.prpd is None or self.waveform is None:
            messagebox.showwarning("Warning", "No data to save.")
            return

        save_dir = filedialog.askdirectory(title="Select Directory to Save Files")
        if not save_dir:
            return

        # Save files using original filenames with '_E' appended before the extension
        prpd_base, prpd_ext = os.path.splitext(self.prpd_filename)
        waveform_base, waveform_ext = os.path.splitext(self.waveform_filename)

        prpd_save_filename = f"{prpd_base}_E{prpd_ext}"
        waveform_save_filename = f"{waveform_base}_E{waveform_ext}"

        prpd_save_path = os.path.join(save_dir, prpd_save_filename)
        waveform_save_path = os.path.join(save_dir, waveform_save_filename)

        try:
            np.savetxt(prpd_save_path, self.prpd, delimiter=',')
            np.savetxt(waveform_save_path, self.waveform, delimiter=',')
            messagebox.showinfo("Success", f"Files saved to {save_dir}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save files:\n{e}")

def main():
    root = tk.Tk()
    app = PDDEA(root)
    root.mainloop()

if __name__ == "__main__":
    main()
