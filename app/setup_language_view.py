"""Reviewable setup editor: compilation never changes physics or runs PDE."""
import json
import tkinter as tk
from tkinter import ttk, messagebox
from .setup_language import emit, compile_setup, analyze_setup


class SetupLanguageView:
    def __init__(self, app):
        self.app = app
        self.window = tk.Toplevel(app.root)
        self.window._aft_title_key = 'setup.title'
        self.window.title(app._tr('setup.title'))
        self.window.geometry('850x600')
        self.text = tk.Text(self.window, wrap='none', undo=True)
        scroll = ttk.Scrollbar(self.window, command=self.text.yview)
        self.text.configure(yscrollcommand=scroll.set)
        scroll.pack(side='right', fill='y')
        buttons = ttk.Frame(self.window); buttons.pack(side='bottom', fill='x')
        self.text.pack(fill='both', expand=True)
        for key, command in [('setup.refresh', self.refresh), ('setup.validate', self.validate),
                             ('setup.apply', self.apply), ('ai.title', self.open_ai)]:
            app._bind_text(ttk.Button(buttons, command=command), key).pack(side='left')
        self.output = tk.StringVar(master=self.window)
        ttk.Label(buttons, textvariable=self.output, wraplength=420).pack(side='left')
        self.refresh()

    def open_ai(self):
        from .ai_chat_view import AIChatView
        AIChatView(self.app, lambda: self.text.get('1.0','end').strip())

    def refresh(self):
        load = self.app.load_workflow
        if load._mesh is None: return
        self.text.delete('1.0', 'end')
        self.text.insert('1.0', emit(load._mesh, load.loads, load.correction))

    def validate(self):
        try:
            mesh = self.app.geometry_workflow.mesh
            if mesh is None: raise ValueError('mesh required')
            loads, correction = compile_setup(self.text.get('1.0','end'), mesh)
            self.output.set(json.dumps(analyze_setup(mesh, loads, correction), ensure_ascii=False))
            return loads, correction
        except (ValueError, TypeError, SyntaxError, KeyError) as exc:
            self.output.set(str(exc)); return None

    def apply(self):
        result = self.validate()
        if result is None: return
        if not messagebox.askyesno(self.app._tr('setup.title'), self.app._tr('setup.confirm'), parent=self.window): return
        loads, correction_faces = result
        workflow = self.app.load_workflow
        workflow.loads = loads
        workflow.applied = loads[-1] if loads else None
        workflow.correction = None
        workflow._refresh_loads()
        if correction_faces is not None:
            workflow.select_faces(correction_faces)
            workflow.balance()
