"""Viewport-only scrolling for long forms; no application/solver state here."""
from __future__ import annotations

import math
import tkinter as tk
from tkinter import ttk


class ScrollablePanel(ttk.Frame):
    """A vertical form viewport with scoped wheel and keyboard bindings.

    Each panel owns a private bindtag on its descendants. No bind_all handler
    steals Matplotlib zoom, and a wheel over a combobox scrolls the form without
    changing that numerical/enum input. Tab focus reveals off-screen fields.
    """

    def __init__(self, master, *, width=320, background="white", **kwargs):
        super().__init__(master, **kwargs)
        self.canvas = tk.Canvas(
            self, width=width, height=1, background=background,
            highlightthickness=0, borderwidth=0, takefocus=True,
            yscrollincrement=18,
        )
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollbar.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        self.content = ttk.Frame(self, style=kwargs.get("style", "TFrame"))
        self._window = self.canvas.create_window(0, 0, window=self.content, anchor="nw")
        self._tag = f"FormScroll:{self}"
        self._wheel_remainder = 0.0
        self._focus_job = None
        self._sequences = {
            "<MouseWheel>": self._wheel,
            "<Button-4>": self._wheel,
            "<Button-5>": self._wheel,
            "<Prior>": lambda _event: self._page(-1),
            "<Next>": lambda _event: self._page(1),
            "<FocusIn>": self._focus,
        }
        for sequence, callback in self._sequences.items():
            self.bind_class(self._tag, sequence, callback)
        self._tag_widget(self.canvas)
        self.canvas.bind("<Home>", lambda _event: self.canvas.yview_moveto(0))
        self.canvas.bind("<End>", lambda _event: self.canvas.yview_moveto(1))
        self.canvas.bind("<Configure>", self._resize)
        self.content.bind("<Configure>", self._layout)

    def _tag_widget(self, widget):
        if self._tag not in widget.bindtags():
            widget.bindtags((self._tag,) + widget.bindtags())
        for child in widget.winfo_children():
            self._tag_widget(child)

    def _layout(self, _event=None):
        self._tag_widget(self.content)
        self.canvas.configure(scrollregion=self.canvas.bbox(self._window))

    def _resize(self, event):
        self.canvas.itemconfigure(self._window, width=event.width)

    def _overflow(self):
        low, high = self.canvas.yview()
        return high - low < 1.0

    def _wheel(self, event):
        if not self._overflow():
            # Still prevent native combobox wheel selection changes.
            return "break"
        if getattr(event, "num", None) in (4, 5):
            increment = -3 if event.num == 4 else 3
        elif self.tk.call("tk", "windowingsystem") == "aqua":
            increment = -event.delta
        else:
            increment = -3 * event.delta / 120
        self._wheel_remainder += increment
        units = math.trunc(self._wheel_remainder)
        self._wheel_remainder -= units
        if units:
            self.canvas.yview_scroll(units, "units")
        return "break"

    def _page(self, direction):
        self.canvas.yview_scroll(direction, "pages")
        return "break"

    def _focus(self, event):
        if event.widget is self.canvas:
            return
        if self._focus_job is not None:
            self.after_cancel(self._focus_job)
        self._focus_job = self.after_idle(lambda: self._reveal(event.widget))

    def _reveal(self, widget):
        self._focus_job = None
        if not widget.winfo_exists() or not self._overflow():
            return
        top = widget.winfo_rooty() - self.content.winfo_rooty()
        bottom = top + widget.winfo_height()
        view_top = self.canvas.canvasy(0)
        view_height = self.canvas.winfo_height()
        target = None
        if top < view_top:
            target = top
        elif bottom > view_top + view_height:
            target = bottom - view_height
        if target is not None:
            self.canvas.yview_moveto(max(0, target) / self.content.winfo_height())

    def destroy(self):
        if self._focus_job is not None:
            self.after_cancel(self._focus_job)
        for sequence in self._sequences:
            self.unbind_class(self._tag, sequence)
        super().destroy()
