from manimlib.imports import *

##states
U = np.array([1.0,  0.0])
D = np.array([0.0,  1.0])
L = np.array([1.0,  1.0]) * 0.5 ** 0.5
R = np.array([1.0, -1.0]) * 0.5 ** 0.5

state_lookup = {
    "u" : U,
    "d" : D,
    "r" : R,
    "l" : L,
    "ud" : (0.5 ** 0.5) * (U + D),
    "lr" : (0.5 ** 0.5) * (L + R)
}

label_lookup = {
    "z" : "Color",
    "x" : "Case"
}

class Funnel(VGroup):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.create_body()

    def create_body(self):
        #corners
        self.ul = ul = UL
        self.dl = dl = DL
        self.ur = ur = RIGHT + UP/3.0
        self.dr = dr = RIGHT + DOWN/3.0
        body = Polygon(ul, dl, dr, ur)
        # r_line = Line(ur, dr)
        # t_line = Line(ur, ul)
        # b_line = Line(dr, dl)
        # l_arc = ArcBetweenPoints(dl, ul)
        self.body = body #VGroup(r_line, b_line, t_line, l_arc)
        self.add(self.body)
        
        # self.add(Dot(self.opening), Dot(self.get_right()))

    def get_opening(self, flipped=False):
        if flipped:
            return self.body.get_right()
        return self.body.get_left()

class Detector(VGroup):

    def __init__(self, axis="z", **kwargs):
        super().__init__(**kwargs)
        self.axis = axis
        self.create_body()
        
    def create_body(self):
        box = Square()
        inp = Funnel().scale(0.3)
        out1 = Funnel().flip().scale(0.3)
        out2 = Funnel().flip().scale(0.3)
        inp.move_to(box.get_left() + inp.get_width() * LEFT / 2.0)
        out1.move_to((box.get_right() + out1.get_width() * RIGHT / 2.0) + box.get_top() / 2.0 * UP * 1.5)
        out2.move_to((box.get_right() + out2.get_width() * RIGHT / 2.0) + box.get_top() / 2.0 * DOWN * 1.5)
        label = TextMobject(label_lookup[self.axis])
        label.set_width(box.get_width() - 0.2)
        label.move_to(self.get_center())
        body = VGroup(inp, out1, out2, box, label)
        self.body = body
        self.inp = inp
        self.out1 = out1
        self.out2 = out2
        # self.inp = inp.get_opening()
        # self.out_top = out1.get_opening(True)
        # self.out_bot = out2.get_opening(True)
        self.add(self.body)

    def get_input_point(self):
        return self.body[0].get_left()

    def get_top_output_point(self):
        return self.body[1].get_right()

    def get_bot_output_point(self):
        return self.body[2].get_right()

    def setup_top_output_destination(self, target):
        if isinstance(target, Detector):
            self.top_output_destination = target.get_input_point()
            self.top_output_destination_type = target
        elif isinstance(target, Collimator):
            self.top_output_destination = target.get_input_point()
            self.top_output_destination_type = target
        elif isinstance(target, Collector):
            self.top_output_destination = target.get_body_center()
            self.top_output_destination_type = target
        else:
            self.top_output_destination = target
            self.top_output_destination_type = None
    
    def setup_bot_output_destination(self, target):
        if isinstance(target, Detector):
            self.bot_output_destination = target.get_input_point()
            self.bot_output_destination_type = target
        elif isinstance(target, Collimator):
            self.bot_output_destination = target.get_input_point()
            self.bot_output_destination_type = target
        elif isinstance(target, Collector):
            self.bot_output_destination = target.get_body_center()
            self.bot_output_destination_type = target
        else:
            self.bot_output_destination = target
            self.bot_output_destination_type = None

        # self.bot_output_destination = target
        # self.bot_output_destination_type = target_type

    def setup_default_output_destinations(self):
        self.top_output_destination = self.get_top_output_point() + RIGHT
        self.top_output_destination_type = None
        self.bot_output_destination = self.get_bot_output_point() + RIGHT
        self.bot_output_destination_type = None

    def measure(self, e):
        new_state = None
        out_point = None
        destination = None
        destination_type = None
        if self.axis == "z":
            prob_d = abs(np.dot(D, e.state))**2
            r = random.random()
            if  r < prob_d:
                new_state = D
                out_point = self.get_top_output_point()
                destination = self.top_output_destination
                destination_type = self.top_output_destination_type
            else:
                new_state = U
                out_point = self.get_bot_output_point()
                destination = self.bot_output_destination
                destination_type = self.bot_output_destination_type

        elif self.axis == "x":
            prob_r = abs(np.dot(R, e.state))**2
            r = random.random()
            if  r < prob_r:
                new_state = R
                out_point = self.get_top_output_point()
                destination = self.top_output_destination
                destination_type = self.top_output_destination_type
            else:
                new_state = L
                out_point = self.get_bot_output_point()
                destination = self.bot_output_destination
                destination_type = self.bot_output_destination_type


        return new_state, out_point, destination, destination_type

class Collimator(VGroup):
    def __init__(self, axis="z", **kwargs):
        super().__init__(**kwargs)
        self.axis = axis
        self.create_body()
        
    def create_body(self):
        box = Square()
        inp = Funnel().scale(0.3)
        out = Funnel().flip().scale(0.3)
        inp.move_to(box.get_left() + inp.get_width() * LEFT / 2.0)
        out.move_to((box.get_right() + out.get_width() * RIGHT / 2.0))
        label = TextMobject("Collimator")
        label.set_width(box.get_width() - 0.2)
        label.move_to(self.get_center())
        body = VGroup(inp, out, box, label)
        self.body = body
        # self.inp = inp.get_opening()
        # self.out_top = out1.get_opening(True)
        # self.out_bot = out2.get_opening(True)
        self.add(self.body)

    def get_input_point(self):
        return self.body[0].get_left()

    def get_output_point(self):
        return self.body[1].get_right()

    def setup_output_destination(self, target):
        if isinstance(target, Detector):
            self.output_destination = target.get_input_point()
            self.output_destination_type = target
        elif isinstance(target, Collimator):
            self.output_destination = target.get_input_point()
            self.output_destination_type = target
        elif isinstance(target, Collector):
            self.output_destination = target.get_body_center()
            self.output_destination_type = target
        else:
            self.output_destination = target
            self.output_destination_type = None
        # self.output_destination = target
        # self.output_destination_type = target_type

    def setup_default_output_destinations(self):
        self.output_destination = self.get_output_point() + RIGHT
        self.output_destination_type = None

    def measure(self, e):
        new_state = None
        out_point = None
        destination = None
        destination_type = None
        if self.axis == "z":
            c1 = np.dot(U, e.state)
            c2 = np.dot(D, e.state)
            new_state = c1 * U + c2 * D
            out_point = self.get_output_point()
            destination = self.output_destination
            destination_type = self.output_destination_type

        elif self.axis == "x":
            c1 = np.dot(R, e.state)
            c2 = np.dot(L, e.state)
            new_state = c1 * R + c2 * L
            out_point = self.get_output_point()
            destination = self.output_destination
            destination_type = self.output_destination_type

        return new_state, out_point, destination, destination_type
        

class Collector(VGroup):
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.create_body()
        self.create_counter()
    
    def create_counter(self):
        self.counter = Integer()
        self.counter.add_updater(lambda m: m.next_to(self.body, RIGHT))
        self.add(self.counter)

    def create_body(self):
        t_line = Line(UL/2, UR/2)
        r_line = Line(UR/2, DR/2)
        b_line = Line(DR/2, DL/2)

        door1 = Line(UL/2, LEFT/2)
        door2 = Line(LEFT/2, DL/2)
        door1.rotate(-(3*PI/2 - PI/6), about_point=UL/2)
        door2.rotate( (3*PI/2 - PI/6), about_point=DL/2)
        body = VGroup(t_line, r_line, b_line)
        self.door = VGroup(door1, door2)
        body.move_to(self.get_center())
        self.body = body
        self.add(self.body)
        self.add(self.door)

    def get_body_center(self):
        return self.body.get_center()

    def update_counter(self):
        self.counter.set_value(self.counter.get_value() + 1)

    def get_number_of_electrons(self):
        return self.counter.get_value()

    def close_door(self):

        #this is why we need alpha!!
        self.door.save_state()
        def closing(mob, alpha):
            mob.restore()
            mob[0].rotate( (3*PI/2 - PI/6) * alpha, about_point=self.body.get_corner(UL))
            mob[1].rotate(-(3*PI/2 - PI/6) * alpha, about_point=self.body.get_corner(DL))
        return UpdateFromAlphaFunc(self.door, closing)

class Electron(VGroup):
    CONFIG = {
        "r" : 0.4
    }

    def __init__(self, state=None, **kwargs):
        super().__init__(**kwargs)
        if state is not None:
            # self.label = state_label
            self.state = state
        else:
            self.state = (U + D + L + R) * 0.5
            # self.label = None
        
        
        self.create_body()
        self.ismoving = False
        self.counted = False
        self.isRandomWalking = False
        self.target_type = None
        self.target_location = self.get_center() + 2 * RIGHT
        self.speed = 3.0
        self.direction = self.target_location - self.get_center()
        self.direction /= np.linalg.norm(self.direction)
        self.timer = 0.0
        self.lag = 0.0
        self.add_updater(lambda m, dt: m.movement(dt))
        

    def create_body(self):
        #U and D are colors
        if abs(np.dot(self.state, U)**2 - 1.0) < EPSILON:
            # body = RoundedRectangle(width=self.r*2, height=self.r*2, corner_radius=self.r/2 - self.r/10, fill_opacity=1.0)
            body = Circle(radius=self.r, fill_opacity=0.7)
            text = TextMobject("??").set_color(YELLOW)
            body.set_color(BLUE)
        elif abs(np.dot(self.state, D)**2 - 1.0) < EPSILON:
            # body = RoundedRectangle(width=self.r*2, height=self.r*2, corner_radius=self.r/2 - self.r/10, fill_opacity=1.0)
            body = Circle(radius=self.r, fill_opacity=0.7)
            text = TextMobject("??").set_color(YELLOW)
            body.set_color(RED)

        #L and R are shapes
        elif abs(np.dot(self.state, R)**2 - 1.0) < EPSILON:
            # body = Rectangle(width=self.r*2, height=self.r*2, fill_opacity=1.0)
            # body.set_color([BLUE, RED])
            # body.set_sheen_direction(RIGHT)
            body = Circle(radius=self.r, fill_opacity=0.0).set_color(ORANGE)
            text = TextMobject("e").set_color(YELLOW)
            
        elif abs(np.dot(self.state, L)**2 - 1.0) < EPSILON:
            # body = Circle(radius=self.r, fill_opacity=1.0)
            # body.set_color([BLUE, RED])
            # body.set_sheen_direction(RIGHT)
            body = Circle(radius=self.r, fill_opacity=0.0).set_color(ORANGE)
            text = TextMobject("E").set_color(YELLOW)
        else:
            # body = RoundedRectangle(width=self.r*2, height=self.r*2, corner_radius=self.r/2 - self.r/10, fill_opacity=1.0)
            # body.set_color([BLUE, RED])
            # body.set_sheen_direction(RIGHT)
            body = Circle(radius=self.r, fill_opacity=0.0).set_color(ORANGE)
            text = TextMobject("??").set_color(YELLOW)
        # self.body = body
        text.move_to(body.get_center())
        self.body = VGroup(body, text)
        self.add(self.body)

    def set_lag(self, lag):
        self.lag = lag

    def trigger_movement(self):
        self.isRandomWalking = False
        self.ismoving = True

    def set_target(self, target, target_type=None):
        self.target_location = target
        self.direction = self.target_location - self.get_center()
        self.direction /= np.linalg.norm(self.direction)
        self.target_type = target_type

    def transform_to_new_state(self, state):
        new_e = Electron(state).move_to(self.get_center()).set_height(self.get_width())
        # self.label = state_label
        self.state = state
        self.become(new_e)

    def get_random_step(self):
        if random.random() < 0.01:
            angle = random.uniform(-PI/2, PI/2)
            mat = np.array(
                [[np.cos(angle), -np.sin(angle), 0.0], 
                [np.sin(angle), np.cos(angle), 0.0], 
                [0., 0., 1.0]]     
                )
            self.direction = np.dot(mat, self.direction)
        pos = self.get_center() 
        if isinstance(self.target_type, Collector) or isinstance(self.target_type, Gun):
            b = self.target_type.body
            r = self.get_width() / 2
            if pos[0] - r < (b.get_left())[0]:
                # self.shift(r * RIGHT)
                self.move_to(np.array([(b.get_left())[0] + r, pos[1], 0]))
                self.direction[0] *= -1
            if pos[0] + r > (b.get_right())[0]:
                # self.shift(r * LEFT)
                self.move_to(np.array([(b.get_right())[0] - r, pos[1], 0]))
                self.direction[0] *= -1
            if pos[1] + r > (b.get_top()) [1]:
                # self.shift(r * DOWN)
                self.move_to(np.array([pos[0], b.get_top()[1] - r, 0]))
                self.direction[1] *= -1
            if pos[1] - r < (b.get_bottom())[1]:
                # self.shift(r * UP)
                self.move_to(np.array([pos[0], b.get_bottom()[1] + r, 0]))
                self.direction[1] *= -1
            
        return (self.speed / 10.0) * self.direction

    def movement(self, dt):
        if self.ismoving and self.timer > self.lag and not self.isRandomWalking:
            c = self.get_center()
            target = self.target_location - c
            d = np.linalg.norm(target)
            if d > self.speed * dt:
                self.shift(self.direction * self.speed * dt)
            elif self.target_type is None: 
                self.ismoving = False
            elif isinstance(self.target_type, Detector):
                new_state, out_point, destination, destination_type = self.target_type.measure(self)
                self.move_to(out_point)
                self.transform_to_new_state(new_state)
                self.set_target(destination, destination_type)
            elif isinstance(self.target_type, Collector):
                self.ismoving = False
                if not self.counted:
                    self.target_type.update_counter()
                    self.target_type.add(self)
                    self.counted = True
                    self.isRandomWalking = True
            elif isinstance(self.target_type, Collimator):
                new_state, out_point, destination, destination_type = self.target_type.measure(self)
                self.move_to(out_point)
                self.transform_to_new_state(new_state)
                self.set_target(destination, destination_type)

            self.timer += dt
        elif self.ismoving and self.timer <= self.lag:
            self.timer += dt

        elif self.isRandomWalking:
            vel = self.get_random_step()
            self.shift(vel * dt)


# class Gun(VGroup):

#     def __init__(self, capacity=5, type_=None, **kwargs):
#         super().__init__(**kwargs)
#         self.capacity = capacity
#         self.type = type_
#         self.create_body()
#         self.create_electrons()

#     def create_body(self):
#         body = Rectangle(width=0.5, height=0.3, fill_opacity=1.0)
#         self.body = body
#         self.body.move_to(self.get_center())
#         self.add(self.body)

#     def get_num_e(self):
#         return self.capacity

#     def create_electrons(self):
#         self.electrons = VGroup()
#         xs = np.linspace(-0.23, 0.23, self.capacity)
#         # ys = np.linspace(-0.13, 0.13, self.capacity)
#         lags = np.linspace(0, 2, self.capacity) + random.uniform(0, 2)
#         for i in range(self.capacity):
#             e = Electron(self.type).scale(self.body.get_width()/5)
#             e.shift(random.choice(xs) * LEFT)
#             e.set_lag(lags[i])
#             self.electrons.add(e)
        
#         # self.electrons.move_to(self.get_center() + )
#         self.add_to_back(self.electrons)

#     def shoot_electrons(self, target):
#         destination = None
#         target_type = None
#         if isinstance(target, Detector):
#             destination = target.get_input_point()
#             target_type = target
#         elif isinstance(target, Collector):
#             destination = target.get_body_center()
#             target_type = target
#         elif isinstance(target, Collimator):
#             destination = target.get_input_point()
#             target_type = target
#         else:
#             destination = target

#         for e in self.electrons:
#             e.set_target(destination, target_type)
#             e.trigger_movement()

class Gun(VGroup):

    def __init__(self, capacity=1, type_=None, **kwargs):
        super().__init__(**kwargs)
        self.capacity = capacity
        self.type = type_
        self.create_body()
        self.create_electrons()

    def create_body(self):
        tline = Line(UL/4, UR/4)
        lline = Line(UL/4, DL/4)
        bline = Line(DL/4, DR/4)
        self.body = VGroup(tline, lline, bline)
        self.add(self.body)
    
    def get_num_e(self):
        return self.capacity

    def create_electrons(self):
        self.electrons = VGroup()
        xs = np.linspace(-0.23, 0.23, self.capacity)
        # ys = np.linspace(-0.13, 0.13, self.capacity)
        lags = np.linspace(0, 2, self.capacity) + random.uniform(0, 2)
        for i in range(self.capacity):
            e = Electron(self.type).scale(self.body.get_width()/8)
            e.shift(random.choice(xs) * LEFT)
            e.set_lag(lags[i])
            e.set_target(self.get_center(), self)
            e.isRandomWalking = True
            self.electrons.add(e)
        
        # self.electrons.move_to(self.get_center() + )
        self.add_to_back(self.electrons)

    def shoot_electrons(self, target):
        destination = None
        target_type = None
        if isinstance(target, Detector):
            destination = target.get_input_point()
            target_type = target
        elif isinstance(target, Collector):
            destination = target.get_body_center()
            target_type = target
        elif isinstance(target, Collimator):
            destination = target.get_input_point()
            target_type = target
        else:
            destination = target

        for e in self.electrons:
            e.set_target(destination, target_type)
            e.trigger_movement()
        

class Test(Scene):
    CONFIG = {
        "random_seed" : None,
        "num_e" : 10
    }
    def construct(self):
        # Det = Detector(axis="x")
        # Det.scale(0.6)
        # self.add(Det.shift(UP + LEFT * 2))
        # Det1 = Detector(axis="z")
        # Det1.scale(0.6)
        # self.add(Det1.shift(Det.get_bot_output_point() + RIGHT * 3))
        # Det.setup_default_output_destinations()
        # Det.setup_top_output_destination(Det1.get_input_point(), Det1)
        # Det.setup_bot_output_destination(Det1.get_input_point(), Det1)
        # Det1.setup_default_output_destinations()
        # self.add(Det1)

        # G = Gun(100, type_=U)
        # self.add(G.move_to(LEFT*5 + UP))
        # G.shoot_electrons(Det.get_input_point(), Det)
        # # col = Collimator("x")
        # # col.setup_default_output_destinations()
        # C = Collector()
        # self.add(C.shift(UP))
        # # self.add(col)
        # # G.shoot_electrons(col.get_input_point(), col)
        # # e = Electron(U)
        # # self.add(e)
        # C = Collector().move_to(UP).scale(0.5)
        # self.add(C, Dot().move_to(C.get_body_center()))
        # self.play(C.close_door())

        G = Gun(self.num_e, U)
        D1 = Detector(axis="z").scale(1)
        D2 = Detector(axis="x").scale(1)
        C1 = Collector().scale(0.4)
        C2 = Collector().scale(0.4)
        C3 = Collector().scale(0.4)

        G.to_edge(LEFT)
        D1.move_to(G.get_center() + RIGHT * 2.5)
        C1.move_to(D1.get_top_output_point() + RIGHT + UP * (C1.get_center() - C1.get_body_center()))
        D2.move_to(D1.get_bot_output_point() + RIGHT * 4)
        C2.move_to(D2.get_top_output_point() + RIGHT*2 + UP * (C2.get_center() - C2.get_body_center()))
        C3.move_to(D2.get_bot_output_point() + RIGHT*2 + UP * (C3.get_center() - C3.get_body_center()))

        D1.setup_top_output_destination(C1)
        D1.setup_bot_output_destination(D2)
        D2.setup_top_output_destination(C2)
        D2.setup_bot_output_destination(C3)
        self.add(G, D1, C1, D2, C2, C3)
        self.wait(5)
        G.shoot_electrons(D1)

        def end_of_expt():
            return C1.get_number_of_electrons() + C2.get_number_of_electrons() + C3.get_number_of_electrons() == G.get_num_e()
        # self.wait(2)
        self.wait_until(end_of_expt)
        self.play(
            C1.close_door(),
            C2.close_door(),
            C3.close_door()
        )
        # self.wait(2)
        # self.remove(G)
        # self.play(
        #     ApplyMethod(C1.scale, 5)
        # )
        # self.play(
        #     FadeOut(VGroup(G, D1, D2, C1)),
        # )
        self.remove(G, D1, D2, C1)
        self.play(
            ApplyMethod(C2.shift, UP*2),
            ApplyMethod(C3.shift, DOWN),
        )
        self.play(
            ApplyMethod(C2.scale, 5),
            ApplyMethod(C3.scale, 5)
        )
    
        self.wait(15)

class IntroElectron(Scene):
    CONFIG = {
        "random_seed": 200
    }
    def construct(self):
        title = Title("e", "l", "e", "c", "t", "r", "o", "n", match_underline_width_to_text=True)
        title.scale(2)
        self.play(
            Write(title)
        )
        p1 = TextMobject("Color").shift(LEFT * 5 + UP * 2).scale(1.5)
        p2 = TextMobject("Case").shift(RIGHT * 5 + UP * 2).scale(1.5)
        self.wait()
        self.play(
            Write(p1)
        )
        self.play(
            Write(p2)
        )
        colors = BulletedList("Blue", "Red", buff=LARGE_BUFF, dot_scale_factor=0)
        colors[0].set_color(BLUE)
        colors[1].set_color(RED)
        colors.to_edge(LEFT)

        cases = BulletedList("Upper", "Lower", buff=LARGE_BUFF, dot_scale_factor=0)
        cases.to_edge(RIGHT)
        self.wait()
        self.play(
            ShowCreation(cases[0]),
            ShowCreation(cases[1])
        )
        self.wait()
        self.play(
            ShowCreation(colors[0]),
            ShowCreation(colors[1])
        )
        # self.add(colors, cases)
        e_U = Electron(U)
        e_D = Electron(D)
        e_R = Electron(R)
        e_L = Electron(L)
        e_U.next_to(colors[0], RIGHT, buff=MED_LARGE_BUFF)
        e_D.next_to(e_U, DOWN, buff=MED_LARGE_BUFF)
        e_R.next_to(cases[1], LEFT, buff=MED_LARGE_BUFF)
        e_L.next_to(cases[0], LEFT, buff=MED_LARGE_BUFF)
        # self.add(e_U, e_D, e_R, e_L)
        self.wait()
        self.play(
            AnimationGroup(
            GrowFromCenter(e_L),
            GrowFromCenter(e_R),
            lag_ratio=1.0,
           )
        )
        self.wait()
        self.play(
            AnimationGroup(
            GrowFromCenter(e_U),
            GrowFromCenter(e_D),
            lag_ratio=1.0,
           )
        )
        self.wait(2)
        self.play(
            AnimationGroup(
                *[FadeOut(x) for x in [p1, p2, e_R, e_L, e_U, e_D, colors, cases]]
            )
        )
        e = Electron()
        self.wait()
        self.play(
            DrawBorderThenFill(e)
        )
        self.wait()
        self.play(
            FadeOut(title)
        )

        self.play(
            ApplyMethod(
                e.shift, LEFT * 4
            )
        )
        self.play(
            ApplyMethod(
                e.scale, 0.7
            )
        )

        colorDet = Detector("z")
        self.play(
            ShowCreation(colorDet)
        )
        inp = colorDet.inp
        out1 = colorDet.out1
        out2 = colorDet.out2

        inp_a = Arrow(inp.get_center() + DL, inp.get_center())
        inp_t = TextMobject("Input").scale(0.7).next_to(inp_a, DL)

        out1_a = Arrow(out1.get_center() + UR, out1.get_center())
        out1_t = TextMobject("Red Output").scale(0.7).next_to(out1_a, UR)

        out2_a = Arrow(out2.get_center() + DR, out2.get_center())
        out2_t = TextMobject("Blue Output").scale(0.7).next_to(out2_a, DR)

        self.play(
            ShowCreation(inp_a), Write(inp_t)
        )
        self.play(
            ShowCreation(out1_a), Write(out1_t),
            ShowCreation(out2_a), Write(out2_t)
        )
        self.wait()
        
        self.play(
            Uncreate(inp_a), Uncreate(inp_t),
            Uncreate(out1_a), Uncreate(out1_t),
            Uncreate(out2_a), Uncreate(out2_t)
        )
        e.save_state()
        colorDet.setup_default_output_destinations()
        e.set_target(colorDet.get_input_point(), colorDet)
        e.trigger_movement()
        self.wait(3)
        self.play(
            ReplacementTransform(colorDet, Detector("x"))
        )
        inp_a = Arrow(inp.get_center() + DL, inp.get_center())
        inp_t = TextMobject("Input").scale(0.7).next_to(inp_a, DL)

        out1_a = Arrow(out1.get_center() + UR, out1.get_center())
        out1_t = TextMobject("CAPS off").scale(0.7).next_to(out1_a, UR)

        out2_a = Arrow(out2.get_center() + DR, out2.get_center())
        out2_t = TextMobject("CAPS on").scale(0.7).next_to(out2_a, DR)
        self.play(
            FadeOut(e)
        )
        self.play(
            ShowCreation(inp_a), Write(inp_t),
            ShowCreation(out1_a), Write(out1_t),
            ShowCreation(out2_a), Write(out2_t)
        )
        self.wait()
        
        self.play(
            Uncreate(inp_a), Uncreate(inp_t),
            Uncreate(out1_a), Uncreate(out1_t),
            Uncreate(out2_a), Uncreate(out2_t)
        )
        colorDet = Detector("x")
        colorDet.setup_default_output_destinations()
        self.play(e.restore)
        e.set_target(colorDet.get_input_point(), colorDet)
        e.trigger_movement()

        self.wait(3)

class Experiment1(Scene):

    CONFIG = {
        "random_seed" : None
    }

    def construct(self):

        G = Gun(1000, (U + D)/2**0.5)
        Det = Detector(axis='z')
        C1 = Collector().scale(0.4)
        C2 = Collector().scale(0.4)

        G.to_edge(LEFT)
        Det.move_to(G.get_center() + RIGHT * 3)
        C1.move_to(Det.get_top_output_point() + RIGHT * 4 + UP * (C1.get_center() - C1.get_body_center()))
        C2.move_to(Det.get_bot_output_point() + RIGHT * 4 + UP * (C2.get_center() - C2.get_body_center()))

        self.add(G, Det, C1, C2)
        self.wait()
        f1 = DashedLine(G.get_right(), Det.get_input_point())
        f2 = DashedLine(Det.get_top_output_point(), C1.get_left())
        f3 = DashedLine(Det.get_bot_output_point(), C2.get_left())

        self.play(
            ShowPassingFlash(f1)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f2),
            ShowPassingFlash(f3)
        )
        self.wait(3)

        Det.setup_bot_output_destination(C2)
        Det.setup_top_output_destination(C1)
        G.shoot_electrons(Det)
        

        def end_of_expt():
            return C1.get_number_of_electrons() + C2.get_number_of_electrons() == G.get_num_e()
        # self.wait(2)
        self.wait_until(end_of_expt)
        self.play(
            C1.close_door(),
            C2.close_door()
        )

        self.remove(G, Det)
        self.play(
            ApplyMethod(C1.shift, UP*2 ),
            ApplyMethod(C2.shift, DOWN ),
        )
        self.play(
            ApplyMethod(C1.scale, 5),
            ApplyMethod(C2.scale, 5)
        )

        e_in_C1 = C1.get_number_of_electrons()
        e_in_C2 = C2.get_number_of_electrons()
        line1 = TextMobject("Out ", "of ", str(e_in_C1 + e_in_C2), " random ", "electrons, ")
        line2 = TextMobject(str(e_in_C2), " are ", "blue ", "and ",  str(e_in_C1), " are ", "red." ).next_to(line1, DOWN)
        result = VGroup(line1, line2).to_edge(LEFT)
                    

        self.play(
            Write(result)
        )
        self.wait(2)

class Experiment2(Scene):

    CONFIG = {
        "random_seed" : None
    }

    def construct(self):

        G = Gun(1000, (U + D)/2**0.5)
        Det1 = Detector(axis='z')
        Det2 = Detector(axis='z')
        C1 = Collector().scale(0.4)
        C2 = Collector().scale(0.4)
        C3 = Collector().scale(0.4)

        G.to_edge(LEFT)
        Det1.move_to(G.get_center() + RIGHT * 3)
        C3.move_to(Det1.get_top_output_point() + RIGHT + UP * (C3.get_center() - C3.get_body_center()))
        Det2.move_to(Det1.get_bot_output_point() + RIGHT * 3)
        C1.move_to(Det2.get_top_output_point() + RIGHT * 3 + UP * (C1.get_center() - C1.get_body_center()))
        C2.move_to(Det2.get_bot_output_point() + RIGHT * 3 + UP * (C2.get_center() - C2.get_body_center()))

        self.add(G, Det1, Det2, C1, C2, C3)
        self.wait()
        f1 = DashedLine(G.get_right(), Det1.get_input_point())
        f2 = DashedLine(Det1.get_top_output_point(), C3.get_left())
        f3 = DashedLine(Det1.get_bot_output_point(), Det2.get_input_point())
        f4 = DashedLine(Det2.get_top_output_point(), C1.get_left())
        f5 = DashedLine(Det2.get_bot_output_point(), C2.get_left())


        self.play(
            ShowPassingFlash(f1)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f2),
            ShowPassingFlash(f3)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f4),
            ShowPassingFlash(f5)
        )
        self.wait(3)

        Det1.setup_bot_output_destination(Det2)
        Det1.setup_top_output_destination(C3)
        Det2.setup_bot_output_destination(C2)
        Det2.setup_top_output_destination(C1)
        G.shoot_electrons(Det1)
        

        def end_of_expt():
            return C1.get_number_of_electrons() + C2.get_number_of_electrons() + C3.get_number_of_electrons() == G.get_num_e()
        # # self.wait(2)
        self.wait_until(end_of_expt)
        self.play(
            C1.close_door(),
            C2.close_door(),
            C3.close_door()
        )

        self.remove(G, Det1, C3, Det2)
        self.play(
            ApplyMethod(C1.shift, UP*2 + LEFT * 0.5),
            ApplyMethod(C2.shift, DOWN + LEFT * 0.5),
        )
        self.play(
            ApplyMethod(C1.scale, 5),
            ApplyMethod(C2.scale, 5)
        )

        e_in_C1 = C1.get_number_of_electrons()
        e_in_C2 = C2.get_number_of_electrons()
        line1 = TextMobject("Out ", "of ", str(e_in_C1 + e_in_C2), " electrons ", "entered in the second detector,")
        line2 = TextMobject(str(e_in_C2), " are ", "blue ", "and ",  str(e_in_C1), " are ", "red." ).next_to(line1, DOWN)
        result = VGroup(line1, line2).to_edge(LEFT)
                    

        self.play(
            Write(result)
        )
        self.wait(2)

class Experiment3(Scene):

    CONFIG = {
        "random_seed" : None
    }

    def construct(self):

        G = Gun(1000, (U + D)/2**0.5)
        Det1 = Detector(axis='z')
        Det2 = Detector(axis='x')
        C1 = Collector().scale(0.4)
        C2 = Collector().scale(0.4)
        C3 = Collector().scale(0.4)

        G.to_edge(LEFT)
        Det1.move_to(G.get_center() + RIGHT * 3)
        C3.move_to(Det1.get_top_output_point() + RIGHT + UP * (C3.get_center() - C3.get_body_center()))
        Det2.move_to(Det1.get_bot_output_point() + RIGHT * 3)
        C1.move_to(Det2.get_top_output_point() + RIGHT * 3 + UP * (C1.get_center() - C1.get_body_center()))
        C2.move_to(Det2.get_bot_output_point() + RIGHT * 3 + UP * (C2.get_center() - C2.get_body_center()))

        self.add(G, Det1, Det2, C1, C2, C3)
        self.wait()
        f1 = DashedLine(G.get_right(), Det1.get_input_point())
        f2 = DashedLine(Det1.get_top_output_point(), C3.get_left())
        f3 = DashedLine(Det1.get_bot_output_point(), Det2.get_input_point())
        f4 = DashedLine(Det2.get_top_output_point(), C1.get_left())
        f5 = DashedLine(Det2.get_bot_output_point(), C2.get_left())


        self.play(
            ShowPassingFlash(f1)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f2),
            ShowPassingFlash(f3)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f4),
            ShowPassingFlash(f5)
        )
        self.wait(3)

        Det1.setup_bot_output_destination(Det2)
        Det1.setup_top_output_destination(C3)
        Det2.setup_bot_output_destination(C2)
        Det2.setup_top_output_destination(C1)
        G.shoot_electrons(Det1)
        

        def end_of_expt():
            return C1.get_number_of_electrons() + C2.get_number_of_electrons() + C3.get_number_of_electrons() == G.get_num_e()
        # # self.wait(2)
        self.wait_until(end_of_expt)
        self.play(
            C1.close_door(),
            C2.close_door(),
            C3.close_door()
        )

        self.remove(G, Det1, C3, Det2)
        self.play(
            ApplyMethod(C1.shift, UP*2 + LEFT * 0.5),
            ApplyMethod(C2.shift, DOWN + LEFT * 0.5),
        )
        self.play(
            ApplyMethod(C1.scale, 5),
            ApplyMethod(C2.scale, 5)
        )

        e_in_C1 = C1.get_number_of_electrons()
        e_in_C2 = C2.get_number_of_electrons()
        line1 = TextMobject("Out ", "of ", str(e_in_C1 + e_in_C2), " electrons ", "entered in the second detector,")
        line2 = TextMobject(str(e_in_C2), " are ", "upper ", "and ",  str(e_in_C1), " are ", "lower." ).next_to(line1, DOWN)
        result = VGroup(line1, line2).to_edge(LEFT)
                    

        self.play(
            Write(result)
        )
        self.wait(2)


class Experiment4(Scene):

    CONFIG = {
        "random_seed" : None
    }

    def construct(self):

        G = Gun(1000, (U + D)/2**0.5)
        Det1 = Detector(axis='z').scale(0.8)
        Det2 = Detector(axis='x').scale(0.8)
        Det3 = Detector(axis='z').scale(0.8)
        C1 = Collector().scale(0.4)
        C2 = Collector().scale(0.4)
        C3 = Collector().scale(0.4)
        C4 = Collector().scale(0.4)

        G.to_edge(LEFT)
        Det1.move_to(G.get_center() + RIGHT * 3)
        C3.move_to(Det1.get_bot_output_point() + RIGHT + UP * (C3.get_center() - C3.get_body_center()))
        Det2.move_to(Det1.get_top_output_point() + RIGHT * 2)
        C4.move_to(Det2.get_bot_output_point() + RIGHT + UP * (C4.get_center() - C4.get_body_center()))
        Det3.move_to(Det2.get_top_output_point() + RIGHT * 2)
        C1.move_to(Det3.get_top_output_point() + RIGHT * 1.5 + UP * (C1.get_center() - C1.get_body_center()))
        C2.move_to(Det3.get_bot_output_point() + RIGHT * 1.5 + UP * (C2.get_center() - C2.get_body_center()))

        self.add(G, Det1, Det2, Det3, C4, C1, C2, C3)
        self.wait()
        f1 = DashedLine(G.get_right(), Det1.get_input_point())
        f2 = DashedLine(Det1.get_bot_output_point(), C3.get_left())
        f3 = DashedLine(Det1.get_top_output_point(), Det2.get_input_point())
        f4 = DashedLine(Det2.get_bot_output_point(), C4.get_left())
        f5 = DashedLine(Det2.get_top_output_point(), Det3.get_input_point())
        f6 = DashedLine(Det3.get_bot_output_point(), C2.get_left())
        f7 = DashedLine(Det3.get_top_output_point(), C1.get_left())


        self.play(
            ShowPassingFlash(f1)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f2),
            ShowPassingFlash(f3)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f4),
            ShowPassingFlash(f5)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f6),
            ShowPassingFlash(f7)
        )
        self.wait(3)

        Det1.setup_top_output_destination(Det2)
        Det1.setup_bot_output_destination(C3)
        Det2.setup_bot_output_destination(C4)
        Det2.setup_top_output_destination(Det3)
        Det3.setup_bot_output_destination(C2)
        Det3.setup_top_output_destination(C1)
        G.shoot_electrons(Det1)
        

        def end_of_expt():
            return C1.get_number_of_electrons() + C2.get_number_of_electrons() + C3.get_number_of_electrons() + C4.get_number_of_electrons() == G.get_num_e()
        # # # self.wait(2)
        self.wait_until(end_of_expt)
        self.play(
            C1.close_door(),
            C2.close_door(),
            C3.close_door(),
            C4.close_door()
        )

        self.remove(G, Det1, C3, Det2, Det3, C4)
        self.play(
            ApplyMethod(C1.shift, UP + LEFT * 0.5),
            ApplyMethod(C2.shift, DOWN * 2 + LEFT * 0.5),
        )
        self.play(
            ApplyMethod(C1.scale, 5),
            ApplyMethod(C2.scale, 5)
        )

        e_in_C1 = C1.get_number_of_electrons()
        e_in_C2 = C2.get_number_of_electrons()
        line1 = TextMobject("Out ", "of ", str(e_in_C1 + e_in_C2), " electrons ", "entered in the third detector,")
        line2 = TextMobject(str(e_in_C2), " are ", "blue ", "and ",  str(e_in_C1), " are ", "red." ).next_to(line1, DOWN)
        result = VGroup(line1, line2).to_edge(LEFT)
                    

        self.play(
            Write(result)
        )
        self.wait(2)

class Experiment5(Scene):

    CONFIG = {
        "random_seed" : None
    }

    def construct(self):

        G = Gun(1000, (U + D)/2**0.5)
        Det1 = Detector(axis='z').scale(0.8)
        Det2 = Collimator(axis='x').scale(0.8)
        Det3 = Detector(axis='z').scale(0.8)
        C1 = Collector().scale(0.4)
        C2 = Collector().scale(0.4)
        C3 = Collector().scale(0.4)
        # C4 = Collector().scale(0.4)

        G.to_edge(LEFT)
        Det1.move_to(G.get_center() + RIGHT * 3)
        C3.move_to(Det1.get_bot_output_point() + RIGHT + UP * (C3.get_center() - C3.get_body_center()))
        Det2.move_to(Det1.get_top_output_point() + RIGHT * 2)
        # C4.move_to(Det2.get_bot_output_point() + RIGHT + UP * (C4.get_center() - C4.get_body_center()))
        Det3.move_to(Det2.get_output_point() + RIGHT * 2)
        C1.move_to(Det3.get_top_output_point() + RIGHT * 1.5 + UP * (C1.get_center() - C1.get_body_center()))
        C2.move_to(Det3.get_bot_output_point() + RIGHT * 1.5 + UP * (C2.get_center() - C2.get_body_center()))

        self.add(G, Det1, Det2, Det3, C1, C2, C3)
        self.wait()
        f1 = DashedLine(G.get_right(), Det1.get_input_point())
        f2 = DashedLine(Det1.get_bot_output_point(), C3.get_left())
        f3 = DashedLine(Det1.get_top_output_point(), Det2.get_input_point())
        f5 = DashedLine(Det2.get_output_point(), Det3.get_input_point())
        f6 = DashedLine(Det3.get_bot_output_point(), C2.get_left())
        f7 = DashedLine(Det3.get_top_output_point(), C1.get_left())


        self.play(
            ShowPassingFlash(f1)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f2),
            ShowPassingFlash(f3)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f5)
        )
        self.wait()
        self.play(
            ShowPassingFlash(f6),
            ShowPassingFlash(f7)
        )
        self.wait(3)

        Det1.setup_top_output_destination(Det2)
        Det1.setup_bot_output_destination(C3)
        Det2.setup_output_destination(Det3)
        Det3.setup_bot_output_destination(C2)
        Det3.setup_top_output_destination(C1)
        G.shoot_electrons(Det1)
        

        def end_of_expt():
            return C1.get_number_of_electrons() + C2.get_number_of_electrons() + C3.get_number_of_electrons() == G.get_num_e()
        # # # self.wait(2)
        self.wait_until(end_of_expt)
        self.play(
            C1.close_door(),
            C2.close_door(),
            C3.close_door()
        )

        self.remove(G, Det1, C3, Det2, Det3)
        self.play(
            ApplyMethod(C1.shift, UP + LEFT * 0.5),
            ApplyMethod(C2.shift, DOWN * 2 + LEFT * 0.5),
        )
        self.play(
            ApplyMethod(C1.scale, 5),
            ApplyMethod(C2.scale, 5)
        )

        e_in_C1 = C1.get_number_of_electrons()
        e_in_C2 = C2.get_number_of_electrons()
        line1 = TextMobject("Out ", "of ", str(e_in_C1 + e_in_C2), " electrons ", "entered in the third detector,")
        line2 = TextMobject(str(e_in_C2), " are ", "blue ", "and ",  str(e_in_C1), " are ", "red." ).next_to(line1, DOWN)
        result = VGroup(line1, line2).to_edge(LEFT)
                    

        self.play(
            Write(result)
        )
        self.wait(2)

class WhatIsCollimator(Scene):

    def construct(self):
        Det = Detector(axis='x').scale(0.5)
        box = Square(side_length=4.0, opacity=0.0)
        inp = Funnel().next_to(box, LEFT, buff=0.0)
        out = Funnel().flip().next_to(box, RIGHT, buff=0.0)

        # self.add(Det, box, inp, out)
        self.play(
            ShowCreation(Det)
        )
        self.wait()
        self.play(
            DrawBorderThenFill(box)
        )
        self.wait()
        self.play(
            ShowCreation(inp),
            ShowCreation(out)
        )
        self.wait()

        f1 = DashedLine(box.get_left(), Det.get_input_point())
        self.play(
            ShowPassingFlash(f1)
        )
        self.wait()
        f2 = DashedLine(Det.get_top_output_point(), box.get_right())
        f3 = DashedLine(Det.get_bot_output_point(), box.get_right())
        self.play(
            ShowPassingFlash(f2),
            ShowPassingFlash(f3)
        )
        self.wait()

class Intro(Scene):
    def construct(self):
        eU = Electron(U)
        eD = Electron(D)
        eR = Electron(R)
        eL = Electron(L)
        e = Electron()

        self.play(
            DrawBorderThenFill(eU)
        )
        self.wait()
        self.play(
            FadeOut(eU)
        )

        self.play(
            DrawBorderThenFill(eD)
        )
        self.wait()
        self.play(
            FadeOut(eD)
        )

        self.play(
            DrawBorderThenFill(eR)
        )
        self.wait()
        self.play(
            FadeOut(eR)
        )

        self.play(
            DrawBorderThenFill(eL)
        )
        self.wait()
        self.play(
            FadeOut(eL)
        )

        self.play(
            DrawBorderThenFill(e)
        )
        self.wait()
        self.play(
            FadeOut(e)
        )
        self.wait()

class IntroMachines(Scene):
    def construct(self):
        e = Electron()
        self.wait()
        self.play(
            DrawBorderThenFill(e)
        )
        self.wait()

        self.play(
            ApplyMethod(
                e.shift, LEFT * 4
            )
        )
        self.play(
            ApplyMethod(
                e.scale, 0.7
            )
        )

        colorDet = Detector("z")
        self.play(
            ShowCreation(colorDet)
        )
        inp = colorDet.inp
        out1 = colorDet.out1
        out2 = colorDet.out2

        inp_a = Arrow(inp.get_center() + DL, inp.get_center())
        inp_t = TextMobject("Input").scale(0.7).next_to(inp_a, DL)

        out1_a = Arrow(out1.get_center() + UR, out1.get_center())
        out1_t = TextMobject("Red Output").scale(0.7).next_to(out1_a, UR)

        out2_a = Arrow(out2.get_center() + DR, out2.get_center())
        out2_t = TextMobject("Blue Output").scale(0.7).next_to(out2_a, DR)

        self.play(
            ShowCreation(inp_a), Write(inp_t)
        )
        self.play(
            ShowCreation(out1_a), Write(out1_t),
            ShowCreation(out2_a), Write(out2_t)
        )
        self.wait()
        
        self.play(
            Uncreate(inp_a), Uncreate(inp_t),
            Uncreate(out1_a), Uncreate(out1_t),
            Uncreate(out2_a), Uncreate(out2_t)
        )
        colorDet.setup_default_output_destinations()
        e.set_target(colorDet.get_input_point(), colorDet)
        e.trigger_movement()
        self.wait(3)
        self.play(
            ReplacementTransform(colorDet, Detector("x"))
        )
        inp_a = Arrow(inp.get_center() + DL, inp.get_center())
        inp_t = TextMobject("Input").scale(0.7).next_to(inp_a, DL)

        out1_a = Arrow(out1.get_center() + UR, out1.get_center())
        out1_t = TextMobject("CAPS off").scale(0.7).next_to(out1_a, UR)

        out2_a = Arrow(out2.get_center() + DR, out2.get_center())
        out2_t = TextMobject("CAPS on").scale(0.7).next_to(out2_a, DR)
        self.play(
            FadeOut(e)
        )
        self.play(
            ShowCreation(inp_a), Write(inp_t),
            ShowCreation(out1_a), Write(out1_t),
            ShowCreation(out2_a), Write(out2_t)
        )
        self.wait()
        
        self.play(
            Uncreate(inp_a), Uncreate(inp_t),
            Uncreate(out1_a), Uncreate(out1_t),
            Uncreate(out2_a), Uncreate(out2_t)
        )
        colorDet = Detector("x")
        colorDet.setup_default_output_destinations()
        e = Electron()
        e.scale(0.7).shift(LEFT * 4)
        self.play(
            FadeIn(e)
        )
        e.set_target(colorDet.get_input_point(), colorDet)
        e.trigger_movement()

        self.wait(3)


# SCENES_IN_ORDER = [
#     Experiment1,
#     Experiment2,
#     Experiment3,
#     Experiment4,
#     Experiment5,
#     WhatIsCollimator
# ]