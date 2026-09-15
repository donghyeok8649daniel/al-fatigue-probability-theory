"""Opt-in, asynchronous advisory chat. No solver or setup mutation hooks."""
import os
import json
import tkinter as tk
from tkinter import ttk
from .aft_ai_controller import ChatController, OpenAIResponsesTransport, Summary, BoundaryError


class AIChatView:
    def __init__(self, app, summary_provider, transport=None):
        self.app, self.summary_provider = app, summary_provider
        self.transport = transport or OpenAIResponsesTransport(lambda: os.environ.get('OPENAI_API_KEY', ''))
        self.is_mock = transport is not None
        self.window = tk.Toplevel(app.root)
        self.window._aft_title_key = 'ai.title'
        self.window.title(app._tr('ai.title'))
        self.window.geometry('900x720')
        self.window.protocol('WM_DELETE_WINDOW', self.close)
        self.controllers = {}
        self.models = {}
        self.current = None
        self.prepared = None
        self.review_state = None
        self.model = tk.StringVar(master=self.window)
        self.reviewed = tk.BooleanVar(master=self.window, value=False)
        self.attach = tk.BooleanVar(master=self.window, value=False)
        self.attach_result = tk.BooleanVar(master=self.window, value=False)
        self.include_history = tk.BooleanVar(master=self.window, value=False)
        self.status = tk.StringVar(master=self.window)
        header = ttk.Frame(self.window); header.pack(fill='x')
        app._bind_text(ttk.Label(header), 'ai.model').pack(side='left')
        ttk.Entry(header, textvariable=self.model, width=28).pack(side='left')
        app._bind_text(ttk.Button(header, command=self.new_session), 'ai.new').pack(side='left')
        self.selector = ttk.Combobox(header, state='readonly', width=15)
        self.selector.pack(side='left')
        self.selector.bind('<<ComboboxSelected>>', self.select_session)
        app._bind_text(ttk.Label(self.window, wraplength=850), 'ai.warning').pack(fill='x')
        for key, variable in [('ai.review',self.reviewed), ('ai.attach',self.attach),
                              ('ai.attach_result',self.attach_result), ('ai.history',self.include_history)]:
            app._bind_text(ttk.Checkbutton(self.window, variable=variable), key).pack(anchor='w')
        self.history = tk.Text(self.window, height=10, wrap='word', state='disabled')
        self.history.pack(fill='both', expand=True)
        self.message = tk.Text(self.window, height=4, wrap='word'); self.message.pack(fill='x')
        self.preview = tk.Text(self.window, height=9, wrap='word', state='disabled')
        self.preview.pack(fill='both', expand=True)
        buttons = ttk.Frame(self.window); buttons.pack(fill='x')
        for key, command in [('ai.preview',self.prepare), ('ai.send',self.send), ('ai.cancel',self.cancel)]:
            app._bind_text(ttk.Button(buttons, command=command), key).pack(side='left')
        ttk.Label(self.window, textvariable=self.status, wraplength=850).pack(fill='x')
        self.timer = self.window.after(100, self.poll)

    @staticmethod
    def set_text(widget, text):
        widget.configure(state='normal'); widget.delete('1.0','end')
        widget.insert('1.0',text); widget.configure(state='disabled')

    def new_session(self):
        try:
            if not self.model.get().strip(): raise BoundaryError('Set an explicit model')
            if len(self.controllers) >= 8: raise BoundaryError('Maximum 8 conversations per window')
            controller = ChatController(self.transport, self.model.get().strip(), max_workers=1)
            sid = controller.new_session()
            name = f'Chat {len(self.controllers)+1}'
            self.controllers[name] = (controller, sid)
            self.models[name] = self.model.get().strip()
            self.selector.configure(values=list(self.controllers))
            self.selector.set(name); self.select_session()
        except BoundaryError as exc: self.status.set(str(exc))

    def select_session(self, _event=None):
        self.current = self.selector.get()
        self.model.set(self.models[self.current])
        self.prepared = None
        self.set_text(self.preview, '')
        self.render_history()

    def render_history(self):
        if self.current not in self.controllers: return
        controller, sid = self.controllers[self.current]
        self.set_text(self.history, '\n\n'.join(f'{m.role}: {m.text}' for m in controller.history(sid)))

    def prepare(self):
        self.prepared = None
        try:
            if self.current not in self.controllers: raise BoundaryError('Create a conversation first')
            if self.model.get().strip() != self.models[self.current]:
                raise BoundaryError('Create a new conversation to change model')
            controller, sid = self.controllers[self.current]
            summaries = (Summary('setup_summary', self.summary_provider()),) if self.attach.get() else ()
            if self.attach_result.get():
                summaries += (Summary('result_summary', self.result_summary()),)
            ids = tuple(m.id for m in controller.history(sid)) if self.include_history.get() else ()
            preview = controller.prepare(sid, self.message.get('1.0','end').strip(), summaries=summaries,
                history_ids=ids, confirmed_non_secret=self.reviewed.get())
            self.prepared = (self.current, preview)
            self.review_state = (self.message.get('1.0','end'), self.attach.get(),
                                 self.include_history.get(), self.model.get(), self.attach_result.get())
            self.set_text(self.preview, preview.payload_json)
            self.status.set(self.app._tr('ai.preview_status'))
        except (BoundaryError, ValueError, tk.TclError) as exc:
            self.set_text(self.preview, ''); self.status.set(str(exc))

    def result_summary(self):
        import numpy as np
        result = self.app.result
        if result is None: raise BoundaryError('No solver result available')
        summary = {'scope':'local_PDE_not_spatial_mechanics', 'crystal_scope':'single_crystal',
                   'experimental_validation':'not_established'}
        for key in ('energy_model_id','time_basis','time_unit','frequency_unit',
                    'probability_resolution_certified','rare_event_floor_scope'):
            if key in result: summary[key] = result[key]
        for key in ('model_time','strain','normal_strain','intrawell_strain','plastic_strain',
                    'cumulative_absorbed_mass','local_rare_event_floor','mass_balance_residual'):
            if key in result:
                array = np.asarray(result[key], dtype=float).reshape(-1)
                if array.size:
                    summary[key] = {'final':float(array[-1]), 'max_abs':float(np.max(np.abs(array)))}
        return json.dumps(summary, ensure_ascii=False, allow_nan=False)

    def send(self):
        try:
            if self.prepared is None or self.prepared[0] != self.current:
                raise BoundaryError('Prepare and review payload first')
            if not self.reviewed.get() or self.review_state != (
                    self.message.get('1.0','end'), self.attach.get(), self.include_history.get(),
                    self.model.get(), self.attach_result.get()):
                raise BoundaryError('Inputs changed; prepare and review again')
            if not self.is_mock and not os.environ.get('OPENAI_API_KEY'):
                raise BoundaryError('OPENAI_API_KEY is not configured')
            controller, _ = self.controllers[self.current]
            controller.send(self.prepared[1], user_clicked_send=True)
            self.prepared = None
            self.status.set(self.app._tr('ai.sending'))
        except BoundaryError as exc: self.status.set(str(exc))

    def cancel(self):
        if self.current in self.controllers:
            controller, sid = self.controllers[self.current]
            controller.cancel(sid); self.prepared = None
            self.status.set(self.app._tr('ai.cancelled'))

    def poll(self):
        for name, (controller, _) in self.controllers.items():
            for event in controller.poll():
                self.status.set(f'{name}: {event.status}')
                if event.status != 'completed': self.status.set(event.text)
                if name == self.current: self.render_history()
        self.timer = self.window.after(100, self.poll)

    def close(self):
        self.window.after_cancel(self.timer)
        for controller, _ in self.controllers.values(): controller.close()
        self.window.destroy()
