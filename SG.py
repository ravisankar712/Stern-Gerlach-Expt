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
    "x" : "Shape"
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
        out1.move_to((box.get_right() + out1.get_width() * RIGHT / 2.0) + box.get_top() / 2.0 * UP)
        out2.move_to((box.get_right() + out2.get_width() * RIGHT / 2.0) + box.get_top() / 2.0 * DOWN)
        label = TextMobject(label_lookup[self.axis])
        label.set_width(box.get_width() - 0.2)
        label.move_to(self.get_center())
        body = VGroup(inp, out1, out2, box, label)
        self.body = body
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

    def setup_top_output_destination(self, target, target_type=None):
        self.top_output_destination = target
        self.top_output_destination_type = target_type
    
    def setup_bot_output_destination(self, target, target_type=None):
        self.bot_output_destination = target
        self.bot_output_destination_type = target_type

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
                new_state = "d"
                out_point = self.get_top_output_point()
                destination = self.top_output_destination
                destination_type = self.top_output_destination_type
            else:
                new_state = "u"
                out_point = self.get_bot_output_point()
                destination = self.bot_output_destination
                destination_type = self.bot_output_destination_type

        elif self.axis == "x":
            prob_r = abs(np.dot(R, e.state))**2
            r = random.random()
            if  r < prob_r:
                new_state = "r"
                out_point = self.get_top_output_point()
                destination = self.top_output_destination
                destination_type = self.top_output_destination_type
            else:
                new_state = "l"
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
        # label = TextMobject(label_lookup[self.axis])
        # label.set_width(box.get_width() - 0.2)
        # label.move_to(self.get_center())
        body = VGroup(inp, out, box)
        self.body = body
        # self.inp = inp.get_opening()
        # self.out_top = out1.get_opening(True)
        # self.out_bot = out2.get_opening(True)
        self.add(self.body)

    def get_input_point(self):
        return self.body[0].get_left()

    def get_output_point(self):
        return self.body[1].get_right()

    def setup_output_destination(self, target, target_type=None):
        self.output_destination = target
        self.output_destination_type = target_type

    def setup_default_output_destinations(self):
        self.output_destination = self.get_output_point() + RIGHT
        self.output_destination_type = None

    def measure(self, e):
        new_state = None
        out_point = None
        destination = None
        destination_type = None
        if self.axis == "z":
            new_state = "ud"
            out_point = self.get_output_point()
            destination = self.output_destination
            destination_type = self.output_destination_type

        elif self.axis == "x":
            new_state = "ud"
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
        body = VGroup(t_line, r_line, b_line)
        body.move_to(self.get_center())
        self.body = body
        self.add(self.body)

    def update_counter(self):
        self.counter.set_value(self.counter.get_value() + 1)


class Electron(VGroup):
    CONFIG = {
        "r" : 0.4
    }

    def __init__(self, state_label=None, **kwargs):
        super().__init__(**kwargs)
        if state_label is not None:
            self.label = state_label
            self.state = state_lookup[state_label]
        else:
            self.state = (U + D + L + R) * 0.5
            self.label = None
        
        
        self.create_body()
        self.ismoving = False
        self.counted = False
        self.isRandomWalking = False
        self.target_type = None
        self.target_location = self.get_center() + 2 * RIGHT
        self.speed = 2.0
        self.direction = self.target_location - self.get_center()
        self.direction /= np.linalg.norm(self.direction)
        self.add_updater(lambda m, dt: m.movement(dt))
        self.timer = 0.0
        self.lag = 0.0

    def create_body(self):
        #U and D are colors
        if self.label == "u":
            body = RoundedRectangle(width=self.r*2, height=self.r*2, corner_radius=self.r/2 - self.r/10, fill_opacity=1.0)
            body.set_color(BLUE)
        elif self.label == "d":
            body = RoundedRectangle(width=self.r*2, height=self.r*2, corner_radius=self.r/2 - self.r/10, fill_opacity=1.0)
            body.set_color(RED)

        #L and R are shapes
        elif self.label == "r":
            body = Rectangle(width=self.r*2, height=self.r*2, fill_opacity=1.0)
            body.set_color([BLUE, RED])
            body.set_sheen_direction(RIGHT)
        elif self.label == "l":
            body = Circle(radius=self.r, fill_opacity=1.0)
            body.set_color([BLUE, RED])
            body.set_sheen_direction(RIGHT)
        else:
            body = RoundedRectangle(width=self.r*2, height=self.r*2, corner_radius=self.r/2 - self.r/10, fill_opacity=1.0)
            body.set_color([BLUE, RED])
            body.set_sheen_direction(RIGHT)
        self.body = body
        self.add(self.body)

    def set_lag(self, lag):
        self.lag = lag

    def trigger_movement(self):
        self.ismoving = True

    def set_target(self, target, target_type=None):
        self.target_location = target
        self.direction = self.target_location - self.get_center()
        self.direction /= np.linalg.norm(self.direction)
        self.target_type = target_type

    def transform_to_new_state(self, state_label):
        new_e = Electron(state_label).move_to(self.get_center()).set_height(self.get_width())
        self.label = state_label
        self.state = state_lookup[state_label]
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
        if isinstance(self.target_type, Collector):
            b = self.target_type.body
            r = self.get_width() / 2
            if pos[0] - r < (b.get_left())[0]:
                self.direction[0] *= -1
            if pos[0] + r > (b.get_right())[0]:
                self.direction[0] *= -1
            if pos[1] + r < (b.get_top()) [1]:
                self.direction[1] *= -1
            if pos[1] - r > (b.get_bottom())[1]:
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


class Gun(VGroup):

    def __init__(self, capacity=5, type_=None, **kwargs):
        super().__init__(**kwargs)
        self.capacity = capacity
        self.type = type_
        self.create_body()
        self.create_electrons()

    def create_body(self):
        body = Rectangle(width=0.5, height=0.3, fill_opacity=1.0)
        self.body = body
        self.body.move_to(self.get_center())
        self.add(self.body)

    def create_electrons(self):
        self.electrons = VGroup()
        xs = np.linspace(-0.25, 0.24, self.capacity)
        lags = np.linspace(0, 1, self.capacity)
        for i in range(self.capacity):
            e = Electron(self.type).scale(self.body.get_width()/5)
            e.shift(xs[i] * LEFT)
            e.set_lag(lags[i])
            self.electrons.add(e)
        
        # self.electrons.move_to(self.get_center() + )
        self.add_to_back(self.electrons)

    def shoot_electrons(self, target, target_type=None):
        for e in self.electrons:
            e.set_target(target, target_type)
            e.trigger_movement()
        

class Test(Scene):
    CONFIG = {
        "random_seed" : None
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

        G = Gun(100, type_="u")
        self.add(G.move_to(LEFT*5 + UP))
        # G.shoot_electrons(Det.get_input_point(), Det)
        col = Collimator("x")
        col.setup_default_output_destinations()
        # C = Collector()
        # self.add(C.shift(UP))
        self.add(col)
        G.shoot_electrons(col.get_input_point(), col)
        self.wait(12)